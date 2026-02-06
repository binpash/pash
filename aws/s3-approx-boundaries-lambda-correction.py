#!/usr/bin/env python3.9
"""
S3 reader that:
- receives approximate chunk boundaries from EC2
- corrects each chunk to line boundaries inside the lambda (no inter-lambda coordination)
"""

import boto3
import os
import sys
import json
import struct
import time

NEWLINE = b"\n"
NEWLINE_INT = 10  # ord("\n")


# ---------------------------------------------------------------------
# Logging / timing helpers
# ---------------------------------------------------------------------

def dprint(debug, msg):
    if debug:
        print(msg, file=sys.stderr)
        sys.stderr.flush()


timing_log = []
timing_start = time.time()


def log_timing(stage, desc="", debug=True):
    elapsed = time.time() - timing_start
    timing_log.append({"stage": stage, "elapsed_ms": int(elapsed * 1000), "desc": desc})
    if debug:
        print(f"[{elapsed:.3f}s][TIMING] {stage}: {desc}", file=sys.stderr)
        sys.stderr.flush()


def print_timing_summary(debug=True):
    if not (debug and timing_log):
        return
    print("\n" + "=" * 70, file=sys.stderr)
    print("[TIMING SUMMARY]", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    prev = 0
    for e in timing_log:
        delta = e["elapsed_ms"] - prev
        print(f"  {e['stage']:20s} @ {e['elapsed_ms']:6d}ms (+{delta:5d}ms) {e['desc']}", file=sys.stderr)
        prev = e["elapsed_ms"]
    print("=" * 70 + "\n", file=sys.stderr)
    sys.stderr.flush()


def _now_ts():
    return f"{time.time():.3f}"


# ---------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------

def parse_keyword_args(args):
    out = {}
    for a in args:
        if "=" in a:
            k, v = a.split("=", 1)
            out[k] = v
    return out


def parse_byte_range(br):
    if not br.startswith("bytes="):
        raise ValueError(f"Invalid byte range format: {br}")
    start_str, end_str = br[6:].split("-", 1)
    return int(start_str), int(end_str)


def parse_chunks(byte_range_str):
    # JSON list mode
    try:
        data = json.loads(byte_range_str)
        if isinstance(data, list) and data:
            return [{
                "start": int(c["start"]),
                "end": int(c["end"]),
                "block_id": int(c["block_id"]),
                "shard_id": int(c.get("shard_id", -1)),
            } for c in data]
    except Exception:
        pass
    # Backwards-compatible single range
    s, e = parse_byte_range(byte_range_str)
    return [{"start": s, "end": e, "block_id": 0, "shard_id": -1}]


# ---------------------------------------------------------------------
# Output header (r_merge compatible)
# ---------------------------------------------------------------------

def write_block_header(fifo, block_id, block_size, is_last):
    # '<qQb7x' => int64 id, uint64 size_t, int8 isLast, 7 pad bytes (total 24 bytes)
    fifo.write(struct.pack("<qQb7x", int(block_id), int(block_size), int(is_last)))


# ---------------------------------------------------------------------
# S3 helper
# ---------------------------------------------------------------------

class S3RangeReader:
    def __init__(self, client, bucket, key):
        self.client = client
        self.bucket = bucket
        self.key = key

    def head_size(self):
        resp = self.client.head_object(Bucket=self.bucket, Key=self.key)
        return int(resp["ContentLength"])

    def get_range(self, start, end):
        if end < start:
            return b""
        resp = self.client.get_object(Bucket=self.bucket, Key=self.key, Range=f"bytes={int(start)}-{int(end)}")
        return resp["Body"].read()

    def get_byte(self, pos):
        if pos < 0:
            return None
        try:
            b = self.get_range(pos, pos)
            return b[0] if b else None
        except Exception:
            return None


# ---------------------------------------------------------------------
# Boundary correction
# ---------------------------------------------------------------------

def resolve_overlap(chunk, window_after_vec, window_size, initial_overlap):
    shard_id = int(chunk.get("shard_id", -1))
    if window_after_vec is not None and shard_id >= 0 and shard_id < len(window_after_vec):
        return window_after_vec[shard_id]
    if window_size is not None:
        return window_size
    return initial_overlap


# ---------------------------------------------------------------------
# Streaming
# ---------------------------------------------------------------------

def stream_chunk_with_correction(
    fifo, s3, chunk, *,
    total_chunks_global,
    file_size,
    window_after_vec,
    window_size,
    initial_overlap,
    write_headers=True,
    debug=False,
):
    start = int(chunk["start"])
    end = int(chunk["end"])
    block_id = int(chunk["block_id"])

    is_last_global = (block_id == total_chunks_global - 1)
    is_last_flag = 1

    dprint(debug, f"[CHUNK {block_id}] bytes={start}-{end}")

    overlap = resolve_overlap(chunk, window_after_vec, window_size, initial_overlap)
    if overlap is None or overlap <= 0:
        overlap = initial_overlap

    read_start = max(0, start - 1)
    if is_last_global:
        read_end = min(file_size - 1, end)
    else:
        read_end = min(file_size - 1, end + overlap)

    buf = bytearray(s3.get_range(read_start, read_end))
    if not buf:
        if write_headers:
            write_block_header(fifo, block_id, 0, is_last_flag)
        fifo.flush()
        return 24 if write_headers else 0

    rel_start = start - read_start
    rel_end = end - read_start
    expansion_count = 0
    max_overlap = 1024 * 1024

    while True:
        # Start alignment
        if block_id == 0:
            adj_rel_start = rel_start
            start_aligned = True
        else:
            if rel_start <= 0:
                adj_rel_start = 0
                start_aligned = True
            elif rel_start - 1 < len(buf) and buf[rel_start - 1] == NEWLINE_INT:
                adj_rel_start = rel_start
                start_aligned = True
            else:
                nl = buf.find(NEWLINE, rel_start)
                if nl != -1:
                    adj_rel_start = nl + 1
                    start_aligned = True
                else:
                    start_aligned = False

        # End alignment
        if is_last_global:
            actual_rel_end = len(buf) - 1
            end_found = True
        else:
            if rel_end >= len(buf):
                end_found = False
            else:
                nl2 = buf.find(NEWLINE, rel_end)
                if nl2 != -1:
                    actual_rel_end = nl2
                    end_found = True
                else:
                    end_found = False

        if start_aligned and end_found:
            break

        if read_end >= file_size - 1:
            if not start_aligned:
                adj_rel_start = len(buf)
            if not end_found:
                actual_rel_end = len(buf) - 1
            break

        # Expand overlap and append tail bytes
        expansion_count += 1
        new_overlap = overlap * 2 if overlap else initial_overlap * 2
        if new_overlap < initial_overlap:
            new_overlap = initial_overlap
        if new_overlap > max_overlap:
            new_overlap = max_overlap

        new_read_end = min(file_size - 1, end + new_overlap)
        if new_read_end <= read_end:
            if not start_aligned:
                adj_rel_start = len(buf)
            if not end_found:
                actual_rel_end = len(buf) - 1
            break

        tail = s3.get_range(read_end + 1, new_read_end)
        if tail:
            buf.extend(tail)
        read_end = new_read_end
        overlap = new_overlap

    if debug:
        dprint(debug, f"[CHUNK {block_id}] overlap={overlap} expansions={expansion_count}")

    if adj_rel_start > actual_rel_end or adj_rel_start >= len(buf) or actual_rel_end < 0:
        if write_headers:
            write_block_header(fifo, block_id, 0, is_last_flag)
        fifo.flush()
        return 24 if write_headers else 0

    data = buf[adj_rel_start:actual_rel_end + 1]
    if write_headers:
        write_block_header(fifo, block_id, len(data), is_last_flag)
    fifo.write(data)
    fifo.flush()
    return (24 + len(data)) if write_headers else len(data)


# ---------------------------------------------------------------------
# Window strategy
# ---------------------------------------------------------------------

def parse_window_strategy(kwargs):
    raw_vec = kwargs.get("window_after_vec")
    window_after_vec = None
    if raw_vec:
        try:
            window_after_vec = [int(x) for x in raw_vec.split(",") if x.strip() != ""]
        except Exception:
            window_after_vec = None

    raw = kwargs.get("window_size", "10240")

    # Prefer shard-specific window vector if provided
    if window_after_vec is not None:
        return window_after_vec, None, 64 * 1024

    # adaptive sentinel without vector -> treat as dynamic initial overlap
    if raw == "adaptive":
        return None, None, 64 * 1024

    # dynamic expansion
    if raw is None or str(raw).lower() in {"none", "null"}:
        return None, None, 64 * 1024

    # fixed
    w = int(raw)
    return None, w, w


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main():
    try:
        log_timing("INIT", "Lambda initialization", debug=True)

        if len(sys.argv) < 4:
            print("Usage: ... <s3_key> <output_fifo> <byte_range> shard=... num_shards=... job_uid=... window_size=...", file=sys.stderr)
            sys.exit(1)

        s3_key, output_fifo, byte_range_str = sys.argv[1], sys.argv[2], sys.argv[3]
        kwargs = parse_keyword_args(sys.argv[4:])
        shard = int(kwargs.get("shard", 0))
        num_shards = int(kwargs.get("num_shards", 1))
        job_uid = kwargs.get("job_uid", "default")

        perf = os.environ.get("PASH_PERF_MODE", "false").lower() == "true"
        debug = kwargs.get("debug", "false").lower() == "true" and not perf

        window_after_vec, window_size, initial_overlap = parse_window_strategy(kwargs)

        bucket = os.environ.get("AWS_BUCKET")
        if not bucket:
            print("Error: AWS_BUCKET environment variable not set", file=sys.stderr)
            sys.exit(1)

        s3_client = boto3.client("s3")
        s3 = S3RangeReader(s3_client, bucket, s3_key)
        file_size = s3.head_size()

        chunks = parse_chunks(byte_range_str)

        log_timing("ARGS_PARSED", f"Parsed {len(chunks)} chunks", debug)

        # total chunks across ALL lambdas (needed for global-last detection)
        cpl = int(kwargs["chunks_per_lambda"]) if "chunks_per_lambda" in kwargs else None
        if cpl is not None:
            total_chunks_global = num_shards * cpl
        else:
            total_chunks_global = max(c["block_id"] for c in chunks) + 1
            print("[WARNING] chunks_per_lambda not provided, using fallback", file=sys.stderr)

        # Parse write_headers flag (default True for backward compatibility)
        write_headers = True
        if "write_headers" in kwargs:
            write_headers_str = kwargs["write_headers"].lower()
            write_headers = (write_headers_str != "false")

        if debug and window_after_vec is not None and cpl is not None and len(window_after_vec) != cpl:
            dprint(debug, f"[WARN] window_after_vec length {len(window_after_vec)} != chunks_per_lambda {cpl}")

        if debug:
            dprint(debug, f"[{_now_ts()}][MAIN {shard}] S3: {bucket}/{s3_key}")
            dprint(debug, f"[{_now_ts()}][MAIN {shard}] File size: {file_size/(1024*1024*1024):.2f} GB")
            dprint(debug, f"[{_now_ts()}][MAIN {shard}] shard: {shard}/{num_shards} UID: {job_uid}")
            dprint(debug, f"[MAIN {shard}] total_chunks_global={total_chunks_global}")

        log_timing("FIFO_OPEN_START", f"Opening FIFO {output_fifo}", debug)
        with open(output_fifo, "wb", buffering=0) as fifo:
            log_timing("FIFO_OPEN_END", "FIFO connected", debug)
            log_timing("STREAM_START", "Streaming chunks", debug)

            total_written = 0
            for chunk in chunks:
                total_written += stream_chunk_with_correction(
                    fifo,
                    s3,
                    chunk,
                    total_chunks_global=total_chunks_global,
                    file_size=file_size,
                    window_after_vec=window_after_vec,
                    window_size=window_size,
                    initial_overlap=initial_overlap,
                    write_headers=write_headers,
                    debug=debug,
                )

            log_timing("STREAM_END", f"Streamed {total_written} bytes", debug)

        log_timing("COMPLETE", "Lambda complete", debug)
        print_timing_summary(debug)

    except Exception as e:
        print(f"[{_now_ts()}][ERROR] Fatal exception in main: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
        sys.exit(1)


if __name__ == "__main__":
    main()
