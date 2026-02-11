#!/usr/bin/env python3.9
"""
S3 chunk reader: pre-aligned boundaries (smart strategy).

Design:
- EC2 pre-computes exact line-aligned chunk ranges.
- Lambda only reads those exact ranges and writes them to the output FIFO.
- No Lambda-to-Lambda coordination.
- `write_headers` controls whether 24-byte r_merge block headers are emitted.
"""

import json
import os
import struct
import sys
import time

import boto3


# Timing infrastructure for performance debugging
TIMING_LOG = []
TIMING_START = time.time()


def log_timing(stage, description="", debug=True):
    elapsed = time.time() - TIMING_START
    TIMING_LOG.append(
        {
            "stage": stage,
            "elapsed_ms": int(elapsed * 1000),
            "timestamp": time.time(),
            "description": description,
        }
    )
    if debug:
        print(f"[{elapsed:.3f}s][TIMING] {stage}: {description}", file=sys.stderr)
        sys.stderr.flush()


def print_timing_summary(debug=True):
    if debug and TIMING_LOG:
        print(f"\n{'=' * 70}", file=sys.stderr)
        print("[TIMING SUMMARY]", file=sys.stderr)
        print(f"{'=' * 70}", file=sys.stderr)
        prev_ms = 0
        for event in TIMING_LOG:
            delta_ms = event["elapsed_ms"] - prev_ms
            print(
                f"  {event['stage']:20s} @ {event['elapsed_ms']:6d}ms "
                f"(+{delta_ms:5d}ms) {event['description']}",
                file=sys.stderr,
            )
            prev_ms = event["elapsed_ms"]
        print(f"{'=' * 70}\n", file=sys.stderr)
        sys.stderr.flush()


def _now_ts():
    return f"{time.time():.3f}"


def parse_keyword_args(args):
    kwargs = {}
    for arg in args:
        if "=" in arg:
            key, value = arg.split("=", 1)
            kwargs[key] = value
    return kwargs


def parse_byte_range(byte_range_str):
    if not byte_range_str.startswith("bytes="):
        raise ValueError(f"Invalid byte range format: {byte_range_str}")
    start_str, end_str = byte_range_str[6:].split("-", 1)
    return int(start_str), int(end_str)


def parse_chunks(byte_range_str):
    """
    Accept either:
    - JSON list of chunk objects, or
    - single `bytes=START-END` range.
    """
    try:
        chunks_data = json.loads(byte_range_str)
        if isinstance(chunks_data, list) and chunks_data:
            chunks = []
            for chunk in chunks_data:
                chunks.append(
                    {
                        "start": int(chunk["start"]),
                        "end": int(chunk["end"]),
                        "block_id": int(chunk.get("block_id", 0)),
                    }
                )
            return chunks
    except (json.JSONDecodeError, ValueError, TypeError, KeyError):
        pass

    start_byte, end_byte = parse_byte_range(byte_range_str)
    return [{"start": start_byte, "end": end_byte, "block_id": 0}]


def write_block_header(fifo_handle, block_id, block_size, is_last):
    # '<qQb7x': int64 id, uint64 block_size, int8 is_last, 7-byte padding
    fifo_handle.write(struct.pack("<qQb7x", int(block_id), int(block_size), int(is_last)))


def stream_chunk(
    fifo_handle,
    s3_client,
    bucket,
    key,
    chunk,
    *,
    write_headers=True,
    debug=False,
):
    start_byte = int(chunk["start"])
    end_byte = int(chunk["end"])
    block_id = int(chunk.get("block_id", 0))

    if end_byte < start_byte:
        payload_size = 0
        if write_headers:
            write_block_header(fifo_handle, block_id, payload_size, 1)
            fifo_handle.flush()
            return 24
        fifo_handle.flush()
        return 0

    payload_size = end_byte - start_byte + 1

    if debug:
        chunk_mb = payload_size / (1024 * 1024)
        print(
            f"[STREAM] block_id={block_id} bytes={start_byte}-{end_byte} "
            f"({chunk_mb:.2f} MB) headers={write_headers}",
            file=sys.stderr,
        )

    if write_headers:
        write_block_header(fifo_handle, block_id, payload_size, 1)

    response = s3_client.get_object(
        Bucket=bucket,
        Key=key,
        Range=f"bytes={start_byte}-{end_byte}",
    )

    # Stream S3 body directly to FIFO in bounded chunks.
    written = 24 if write_headers else 0
    for body_chunk in response["Body"].iter_chunks(chunk_size=8 * 1024 * 1024):
        if body_chunk:
            fifo_handle.write(body_chunk)
            written += len(body_chunk)

    fifo_handle.flush()
    return written


def main():
    try:
        log_timing("INIT", "Lambda initialization", debug=True)

        if len(sys.argv) < 4:
            print(
                "Usage: python3.9 aws/s3-chunk-reader-smart-prealigned.py "
                "<s3_key> <output_fifo> <byte_range> shard=<N> num_shards=<N> "
                "job_uid=<UID> [write_headers=true|false]",
                file=sys.stderr,
            )
            sys.exit(1)

        s3_key = sys.argv[1]
        output_fifo = sys.argv[2]
        byte_range_str = sys.argv[3]

        kwargs = parse_keyword_args(sys.argv[4:])
        shard_index = int(kwargs.get("shard", 0))
        num_shards = int(kwargs.get("num_shards", 1))
        job_uid = kwargs.get("job_uid", "default")

        perf_mode = os.environ.get("PASH_PERF_MODE", "false").lower() == "true"
        debug = kwargs.get("debug", "false").lower() == "true" and not perf_mode

        write_headers = kwargs.get("write_headers", "true").lower() != "false"

        s3_bucket = os.environ.get("AWS_BUCKET")
        if not s3_bucket:
            print("Error: AWS_BUCKET environment variable not set", file=sys.stderr)
            sys.exit(1)

        chunks = parse_chunks(byte_range_str)

        if debug:
            chunk_ids = [chunk["block_id"] for chunk in chunks]
            print(f"[{_now_ts()}][MAIN {shard_index}] strategy=s3_smart_prealigned", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] shard: {shard_index}/{num_shards}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] UID: {job_uid}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] chunks: {chunk_ids}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] write_headers: {write_headers}", file=sys.stderr)
            sys.stderr.flush()

        s3_client = boto3.client("s3")

        log_timing("FIFO_OPEN_START", f"Opening FIFO {output_fifo}", debug)
        with open(output_fifo, "wb", buffering=0) as fifo:
            log_timing("FIFO_OPEN_END", "FIFO connected", debug)
            log_timing("STREAM_START", f"Streaming {len(chunks)} chunk(s)", debug)

            total_written = 0
            for chunk in chunks:
                total_written += stream_chunk(
                    fifo,
                    s3_client,
                    s3_bucket,
                    s3_key,
                    chunk,
                    write_headers=write_headers,
                    debug=debug,
                )

            log_timing("STREAM_END", f"Streamed {total_written} bytes", debug)

        log_timing("COMPLETE", "Lambda complete", debug)
        print_timing_summary(debug)

    except Exception as exc:
        print(f"[{_now_ts()}][ERROR] Fatal exception in main: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
        sys.exit(1)


if __name__ == "__main__":
    main()
