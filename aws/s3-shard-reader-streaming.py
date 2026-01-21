#!/usr/bin/env python3.9
"""
S3 Shard Reader with TRUE Streaming and Lambda-to-Lambda Communication

Key difference from s3-shard-reader.py (smart boundaries mode):
- Streams data to FIFO immediately as it downloads
- Downstream shell commands start processing right away
- Only sequential dependency is for the LAST LINE (tail byte request)
- All other lines are processed in parallel across all lambdas

Architecture:
  Lambda 0: [====== stream all lines ======][wait for Lambda 1 to request tail]
  Lambda 1: [====== stream all lines ======][request tail from 0][wait for 2]
  Lambda 2: [====== stream all lines ======][request tail from 1][wait for 3]
  ...

Usage:
    python3.9 s3-shard-reader-streaming.py <s3_key> <output_fifo> <byte_range> \
        shard=<N> num_shards=<N> job_uid=<UID>
"""

import boto3
import os
import sys
import subprocess
import threading
import time
import shutil
from botocore.exceptions import ClientError


# Timing infrastructure for performance debugging
timing_log = []
timing_start = time.time()

def log_timing(stage, description="", debug=True):
    """Log timing event with millisecond precision"""
    elapsed = time.time() - timing_start
    event = {
        "stage": stage,
        "elapsed_ms": int(elapsed * 1000),
        "timestamp": time.time(),
        "description": description
    }
    timing_log.append(event)
    if debug:
        print(f"[{elapsed:.3f}s][TIMING] {stage}: {description}", file=sys.stderr)
        sys.stderr.flush()

def print_timing_summary(debug=True):
    """Print timing summary at end"""
    if debug and timing_log:
        print(f"\n{'='*70}", file=sys.stderr)
        print(f"[TIMING SUMMARY]", file=sys.stderr)
        print(f"{'='*70}", file=sys.stderr)
        prev_ms = 0
        for event in timing_log:
            delta_ms = event['elapsed_ms'] - prev_ms
            print(f"  {event['stage']:20s} @ {event['elapsed_ms']:6d}ms (+{delta_ms:5d}ms) {event['description']}", file=sys.stderr)
            prev_ms = event['elapsed_ms']
        print(f"{'='*70}\n", file=sys.stderr)
        sys.stderr.flush()


def _now_ts():
    return f"{time.time():.3f}"


def _get_env_timeout(name, default_seconds):
    raw = os.environ.get(name)
    if not raw:
        return default_seconds
    try:
        value = float(raw)
    except ValueError:
        return default_seconds
    return max(0.0, value)


RENDEZVOUS_PREFIX = os.environ.get("PASH_S3_RENDEZVOUS_PREFIX", "rendezvous")

# Configurable chunk sizes
CHUNK_SIZE = int(os.environ.get("PASH_S3_CHUNK_SIZE", str(8*1024*1024)))  # 8MB default
RECV_READ_SIZE = int(os.environ.get("PASH_S3_RECV_BUFFER", str(64*1024)))  # 64KB default


def _rendezvous_key(job_uid, requester_index, target_index, kind):
    return f"{RENDEZVOUS_PREFIX}/{job_uid}/{kind}_{requester_index}_to_{target_index}"


def _s3_marker_exists(s3, bucket, key):
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        code = e.response.get("Error", {}).get("Code", "")
        if code in ("404", "NoSuchKey", "NotFound"):
            return False
        raise


def _s3_put_marker(s3, bucket, key, body=b"1"):
    s3.put_object(Bucket=bucket, Key=key, Body=body)


def _s3_delete_marker(s3, bucket, key):
    try:
        s3.delete_object(Bucket=bucket, Key=key)
    except ClientError:
        pass


def _wait_for_marker(s3, bucket, key, timeout_secs, poll_secs, debug=False, log_prefix=""):
    deadline = time.time() + timeout_secs
    while True:
        if _s3_marker_exists(s3, bucket, key):
            if debug and log_prefix:
                print(f"[{_now_ts()}]{log_prefix} Found marker: {key}", file=sys.stderr)
            return True
        if time.time() >= deadline:
            if debug and log_prefix:
                print(f"[{_now_ts()}]{log_prefix} Timeout waiting for marker: {key}", file=sys.stderr)
            return False
        time.sleep(poll_secs)


# Dynamic pashlib path detection
PASHLIB_PATH = os.environ.get('PASHLIB_PATH')
if not PASHLIB_PATH:
    if os.path.exists('/opt/pashlib'):
        PASHLIB_PATH = '/opt/pashlib'
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        PASHLIB_PATH = os.path.join(os.path.dirname(script_dir), 'runtime', 'pashlib')

if not os.path.exists(PASHLIB_PATH):
    raise FileNotFoundError(f"Pashlib not found at {PASHLIB_PATH}")


def parse_byte_range(byte_range_str):
    if not byte_range_str.startswith('bytes='):
        raise ValueError(f"Invalid byte range format: {byte_range_str}")
    range_part = byte_range_str[6:]
    start_str, end_str = range_part.split('-')
    return int(start_str), int(end_str)


def parse_keyword_args(args):
    kwargs = {}
    for arg in args:
        if '=' in arg:
            key, value = arg.split('=', 1)
            kwargs[key] = value
    return kwargs


def stream_s3_to_fifo(fifo_handle, bucket, key, start_byte, end_byte, debug=False):
    """
    Stream S3 byte range DIRECTLY to the FIFO using shutil.copyfileobj.

    Simple approach: no line tracking, just direct streaming.
    Returns number of bytes streamed.
    """
    s3 = boto3.client('s3')

    if debug:
        print(f"[{_now_ts()}][STREAM] Direct streaming: bytes={start_byte}-{end_byte}", file=sys.stderr)
        sys.stderr.flush()

    # Get S3 object with byte range
    response = s3.get_object(
        Bucket=bucket,
        Key=key,
        Range=f'bytes={start_byte}-{end_byte}'
    )

    # Directly pipe S3 stream to FIFO
    # copyfileobj reads chunks and writes them - simple and efficient
    shutil.copyfileobj(response['Body'], fifo_handle, length=CHUNK_SIZE)

    bytes_streamed = end_byte - start_byte + 1

    if debug:
        print(f"[{_now_ts()}][STREAM] Streamed {bytes_streamed} bytes", file=sys.stderr)
        sys.stderr.flush()

    return bytes_streamed


def skip_first_partial_line(bucket, key, start_byte, end_byte, debug=False):
    """
    For shard > 0 with approximate boundaries, skip bytes from start_byte
    until (and including) the first newline.

    Returns: new_start_byte (position after first newline)
    """
    s3 = boto3.client('s3')
    bytes_read = 0
    read_buffer_size = min(CHUNK_SIZE, end_byte - start_byte + 1)

    if debug:
        print(f"[{_now_ts()}][SKIP] Looking for first newline starting at byte {start_byte}", file=sys.stderr)
        sys.stderr.flush()

    while start_byte + bytes_read <= end_byte:
        # Read next chunk
        chunk_end = min(start_byte + bytes_read + read_buffer_size - 1, end_byte)
        response = s3.get_object(
            Bucket=bucket,
            Key=key,
            Range=f'bytes={start_byte + bytes_read}-{chunk_end}'
        )

        chunk = response['Body'].read()
        if not chunk:
            # Reached end without finding newline
            if debug:
                print(f"[{_now_ts()}][SKIP] No newline found, returning end_byte", file=sys.stderr)
                sys.stderr.flush()
            return end_byte

        # Look for newline in chunk
        try:
            newline_pos = chunk.index(b'\n')
            new_start = start_byte + bytes_read + newline_pos + 1
            if debug:
                print(f"[{_now_ts()}][SKIP] Found newline at byte {new_start - 1}, skipped {new_start - start_byte} bytes", file=sys.stderr)
                sys.stderr.flush()
            return new_start
        except ValueError:
            # No newline in this chunk, continue reading
            bytes_read += len(chunk)

    # Reached end_byte without finding newline
    if debug:
        print(f"[{_now_ts()}][SKIP] No newline found in range, returning end_byte", file=sys.stderr)
        sys.stderr.flush()
    return end_byte


def start_tail_server_thread(shard_index, job_uid, s3_bucket, bucket, key, start_byte, end_byte, debug=False):
    """
    Background thread that serves tail bytes to requesting lambdas.

    IMPORTANT: This runs IN PARALLEL with streaming, not after.
    The server sends data from the START of this shard (for the previous lambda's tail).
    """

    def server():
        s3 = boto3.client('s3')
        poll_secs = max(0.05, _get_env_timeout("PASH_S3_TAIL_DISCOVERY_POLL_SECS", 0.2))
        served = set()

        while True:
            for requester_index in range(shard_index):
                if requester_index in served:
                    continue

                req_key = _rendezvous_key(job_uid, requester_index, shard_index, "req")
                if _s3_marker_exists(s3, s3_bucket, req_key):
                    if debug:
                        print(f"[{_now_ts()}][SERVER {shard_index}] Request from {requester_index}", file=sys.stderr)

                    # Serve this requester
                    _serve_tail_bytes(shard_index, requester_index, job_uid, s3_bucket,
                                     bucket, key, start_byte, end_byte, debug, req_key)
                    served.add(requester_index)

            if len(served) == shard_index:
                break
            time.sleep(poll_secs)

    thread = threading.Thread(target=server, daemon=True)
    thread.start()
    return thread


def _serve_tail_bytes(shard_index, requester_index, job_uid, s3_bucket, bucket, key, start_byte, end_byte, debug, req_key):
    """Serve tail bytes from START of this shard to requester."""
    try:
        comm_uid = f"{job_uid}-split-{requester_index}-to-{shard_index}"
        send_fifo = f"/tmp/pash_send_{shard_index}_from_{requester_index}_{os.getpid()}"
        ready_key = _rendezvous_key(job_uid, requester_index, shard_index, "ready")

        try:
            os.unlink(send_fifo)
        except FileNotFoundError:
            pass
        os.mkfifo(send_fifo)

        # Start pashlib send
        pashlib_cmd = [PASHLIB_PATH, f"send*{comm_uid}*0*1*{send_fifo}"]
        proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)

        # Signal ready
        s3 = boto3.client('s3')
        _s3_put_marker(s3, s3_bucket, ready_key, b"ready")

        if debug:
            print(f"[{_now_ts()}][SERVER {shard_index}] Ready for {requester_index}, streaming from S3", file=sys.stderr)

        # Stream bytes from START of this shard until newline
        fifo_fd = os.open(send_fifo, os.O_WRONLY)
        with os.fdopen(fifo_fd, 'wb', buffering=0) as f:
            bytes_sent = 0

            # Read from S3 starting at this shard's start_byte
            response = boto3.client('s3').get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={start_byte}-{end_byte}'
            )

            for chunk in response['Body'].iter_chunks(chunk_size=CHUNK_SIZE):
                try:
                    newline_pos = chunk.index(b'\n')
                    # Found newline - send up to and including it, then stop
                    f.write(chunk[:newline_pos + 1])
                    bytes_sent += newline_pos + 1
                    if debug:
                        print(f"[{_now_ts()}][SERVER {shard_index}] Sent {bytes_sent} bytes to {requester_index} (found newline)", file=sys.stderr)
                    break
                except ValueError:
                    # No newline - send entire chunk
                    f.write(chunk)
                    bytes_sent += len(chunk)

        proc.wait()
        os.unlink(send_fifo)

        if debug:
            print(f"[{_now_ts()}][SERVER {shard_index}] Complete serving {requester_index}", file=sys.stderr)

    except Exception as e:
        print(f"[{_now_ts()}][SERVER {shard_index}] Error serving {requester_index}: {e}", file=sys.stderr)
    finally:
        _s3_delete_marker(s3, s3_bucket, req_key)
        _s3_delete_marker(s3, s3_bucket, ready_key)


def request_tail_bytes(shard_index, num_shards, job_uid, s3_bucket, debug=False):
    """
    Request tail bytes from next shard(s) to complete the last line.

    This is called AFTER streaming is complete, only for the last line.
    """
    accumulated_bytes = bytearray()
    s3 = boto3.client('s3')
    discovery_timeout = _get_env_timeout("PASH_S3_TAIL_DISCOVERY_TIMEOUT_SECS", 30.0)
    poll_secs = max(0.05, _get_env_timeout("PASH_S3_TAIL_DISCOVERY_POLL_SECS", 0.2))

    for target_shard in range(shard_index + 1, num_shards):
        comm_uid = f"{job_uid}-split-{shard_index}-to-{target_shard}"
        req_key = _rendezvous_key(job_uid, shard_index, target_shard, "req")
        ready_key = _rendezvous_key(job_uid, shard_index, target_shard, "ready")

        _s3_put_marker(s3, s3_bucket, req_key, b"req")
        if debug:
            print(f"[{_now_ts()}][CLIENT {shard_index}] Requested tail from shard {target_shard}", file=sys.stderr)

        if not _wait_for_marker(s3, s3_bucket, ready_key, discovery_timeout, poll_secs, debug, f"[CLIENT {shard_index}]"):
            if debug:
                print(f"[{_now_ts()}][CLIENT {shard_index}] No reply from {target_shard}, stopping", file=sys.stderr)
            _s3_delete_marker(s3, s3_bucket, req_key)
            break

        recv_fifo = f"/tmp/pash_recv_{shard_index}_to_{target_shard}_{os.getpid()}"
        try:
            os.unlink(recv_fifo)
        except FileNotFoundError:
            pass
        os.mkfifo(recv_fifo)

        try:
            pashlib_cmd = [PASHLIB_PATH, f"recv*{comm_uid}*1*0*{recv_fifo}"]
            proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)

            found_newline = False
            recv_fd = os.open(recv_fifo, os.O_RDONLY)
            try:
                while True:
                    chunk = os.read(recv_fd, RECV_READ_SIZE)
                    if not chunk:
                        break

                    try:
                        newline_pos = chunk.index(b'\n')
                        accumulated_bytes.extend(chunk[:newline_pos + 1])
                        found_newline = True
                        break
                    except ValueError:
                        accumulated_bytes.extend(chunk)
            finally:
                os.close(recv_fd)

            proc.wait()

            if found_newline:
                if debug:
                    print(f"[{_now_ts()}][CLIENT {shard_index}] Got {len(accumulated_bytes)} tail bytes", file=sys.stderr)
                break

        finally:
            try:
                os.unlink(recv_fifo)
            except FileNotFoundError:
                pass
            _s3_delete_marker(s3, s3_bucket, req_key)
            _s3_delete_marker(s3, s3_bucket, ready_key)

    return bytes(accumulated_bytes)


def main():
    try:
        # Log initialization
        log_timing("INIT", "Lambda initialization", debug=True)

        if len(sys.argv) < 4:
            print("Usage: python3.9 s3-shard-reader-streaming.py <s3_key> <output_fifo> <byte_range> shard=<N> num_shards=<N> job_uid=<UID>", file=sys.stderr)
            sys.exit(1)

        s3_key = sys.argv[1]
        output_fifo = sys.argv[2]
        byte_range_str = sys.argv[3]

        kwargs = parse_keyword_args(sys.argv[4:])
        shard_index = int(kwargs.get('shard', 0))
        num_shards = int(kwargs.get('num_shards', 1))
        job_uid = kwargs.get('job_uid', 'default')
        skip_first_line = kwargs.get('skip_first_line', 'true').lower() == 'true'

        # PERF_MODE disables debug output for clean performance tests
        PERF_MODE = os.environ.get('PASH_PERF_MODE', 'false').lower() == 'true'
        debug = kwargs.get('debug', 'false').lower() == 'true' and not PERF_MODE

        s3_bucket = os.environ.get('AWS_BUCKET')
        if not s3_bucket:
            print("Error: AWS_BUCKET environment variable not set", file=sys.stderr)
            sys.exit(1)

        start_byte, end_byte = parse_byte_range(byte_range_str)
        log_timing("ARGS_PARSED", f"byte_range={start_byte}-{end_byte}", debug)

        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Starting (SIMPLIFIED STREAMING MODE)", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] Range: {byte_range_str}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] shard: {shard_index}/{num_shards}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] skip_first_line: {skip_first_line and shard_index > 0}", file=sys.stderr)
            sys.stderr.flush()

        # For shard > 0: skip first partial line
        actual_start_byte = start_byte
        if skip_first_line and shard_index > 0:
            log_timing("SKIP_START", "Looking for first newline", debug)
            actual_start_byte = skip_first_partial_line(s3_bucket, s3_key, start_byte, end_byte, debug)
            log_timing("SKIP_END", f"Adjusted start: {start_byte} -> {actual_start_byte}", debug)

        # Open FIFO FIRST (enables downstream to start consuming immediately)
        log_timing("FIFO_OPEN_START", f"Opening FIFO {output_fifo}", debug)
        with open(output_fifo, 'wb', buffering=0) as fifo:
            log_timing("FIFO_OPEN_END", "FIFO connected to downstream", debug)
            if debug:
                print(f"[{_now_ts()}][MAIN {shard_index}] FIFO connected, starting parallel operations", file=sys.stderr)
                sys.stderr.flush()

            # Start tail server in PARALLEL (for previous lambdas requesting our data)
            if shard_index > 0:
                if debug:
                    print(f"[{_now_ts()}][MAIN {shard_index}] Starting tail server for {shard_index} potential requesters", file=sys.stderr)
                start_tail_server_thread(shard_index, job_uid, s3_bucket, s3_bucket, s3_key,
                                        start_byte, end_byte, debug)

            # Simple direct streaming
            log_timing("STREAM_START", "Starting S3→FIFO stream", debug)
            bytes_streamed = stream_s3_to_fifo(fifo, s3_bucket, s3_key, actual_start_byte, end_byte, debug)
            log_timing("STREAM_END", f"Streamed {bytes_streamed} bytes", debug)

            if debug:
                print(f"[{_now_ts()}][MAIN {shard_index}] S3 streaming complete: {bytes_streamed} bytes", file=sys.stderr)
                sys.stderr.flush()

            # Request tail bytes from NEXT shard (only for last line!)
            if shard_index < num_shards - 1:
                if debug:
                    print(f"[{_now_ts()}][MAIN {shard_index}] Requesting tail bytes for last line", file=sys.stderr)

                log_timing("TAIL_REQUEST_START", f"Requesting tail from shard {shard_index + 1}", debug)
                tail_bytes = request_tail_bytes(shard_index, num_shards, job_uid, s3_bucket, debug)
                log_timing("TAIL_REQUEST_END", f"Received {len(tail_bytes)} tail bytes", debug)

                if tail_bytes:
                    fifo.write(tail_bytes)
                    fifo.flush()

                    if debug:
                        print(f"[{_now_ts()}][MAIN {shard_index}] Wrote {len(tail_bytes)} tail bytes", file=sys.stderr)

        log_timing("FIFO_CLOSE", "Closing FIFO", debug)
        log_timing("COMPLETE", "Lambda complete", debug)
        print_timing_summary(debug)

        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Complete", file=sys.stderr)
            sys.stderr.flush()

    except Exception as e:
        print(f"[{_now_ts()}][ERROR] Fatal exception: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
