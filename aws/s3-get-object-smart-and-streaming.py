#!/usr/bin/env python3.9
"""
Phase 2: S3 split Reader with Pashlib Inter-Process Communication and Multi-Chunk Headers

Reads one logical split from S3, handling newline boundaries by communicating
with adjacent processes via pashlib FIFOs.

Ports the core algorithm from client/dfs_split_reader.go:
- Lines 96-101: skip_until_newline()
- Lines 62-82: read_until_newline() via pashlib
- Lines 111-146: main orchestration logic

Usage:
    python3.9 s3-shard-reader.py <s3_key> <output_fifo> <byte_range> \
        split=<N> num_shards=<N> job_uid=<UID>

Example:
    python3.9 s3-shard-reader.py "oneliners/inputs/1M.txt" /tmp/fifo21 \
        bytes=0-524287 split=0 num_shards=4 job_uid=test-123
"""

import boto3
import os
import sys
import subprocess
import threading
import time
import shutil
import json
from botocore.exceptions import ClientError
from io import BytesIO


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


# Configurable chunk sizes for performance tuning
CHUNK_SIZE = int(os.environ.get("PASH_S3_CHUNK_SIZE", str(8*1024*1024)))  # Default 8MB
FIFO_BUFFER_SIZE = int(os.environ.get("PASH_S3_FIFO_BUFFER", str(8*1024*1024)))  # Default 8MB
RECV_READ_SIZE = int(os.environ.get("PASH_S3_RECV_BUFFER", str(4*1024)))  # Default 4KB



# Dynamic pashlib path detection
# Works in: Lambda (/opt/pashlib), Local (runtime/pashlib), Custom (PASHLIB_PATH env)
PASHLIB_PATH = os.environ.get('PASHLIB_PATH')
if not PASHLIB_PATH:
    if os.path.exists('/opt/pashlib'):
        PASHLIB_PATH = '/opt/pashlib'  # Lambda environment
    else:
        # Local testing - relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        PASHLIB_PATH = os.path.join(os.path.dirname(script_dir), 'runtime', 'pashlib')

# Verify pashlib exists
if not os.path.exists(PASHLIB_PATH):
    raise FileNotFoundError(f"Pashlib not found at {PASHLIB_PATH}")


def parse_byte_range(byte_range_str):
    """
    Parse byte range string "bytes=0-524287" into (start, end).
    Returns (start_byte, end_byte) inclusive.
    """
    if not byte_range_str.startswith('bytes='):
        raise ValueError(f"Invalid byte range format: {byte_range_str}")

    range_part = byte_range_str[6:]  # Remove 'bytes='
    start_str, end_str = range_part.split('-')
    return int(start_str), int(end_str)


def parse_keyword_args(args):
    """
    Parse keyword arguments like ['split=0', 'num_shards=4', 'job_uid=test']
    Returns dict {'split': '0', 'num_shards': '4', 'job_uid': 'test'}
    """
    kwargs = {}
    for arg in args:
        if '=' in arg:
            key, value = arg.split('=', 1)
            kwargs[key] = value
    return kwargs




def stream_s3_to_fifo_handle(fifo_handle, bucket, key, chunks, debug=False):
    """
    Stream multiple S3 byte ranges to FIFO with r_merge-compatible headers.

    Args:
        fifo_handle: File handle opened in binary write mode
        bucket: S3 bucket name
        key: S3 object key
        chunks: List of dicts with {start, end, block_id, skip_first_line}
        debug: Enable debug output

    Returns:
        Total bytes written (headers + data)
    """
    s3 = boto3.client('s3')
    total_bytes_written = 0

    if debug:
        print(f"[STREAM] Processing {len(chunks)} chunks with headers", file=sys.stderr)

    for chunk_idx, chunk in enumerate(chunks):
        start_byte = chunk['start']
        end_byte = chunk['end']
        block_id = chunk['block_id']
        is_last = 1 #if chunk_idx == len(chunks) - 1 else 0  # Last chunk for THIS lambda

        # Download chunk from S3
        if debug:
            chunk_mb = (end_byte - start_byte + 1) / (1024 * 1024)
            print(f"[STREAM] Chunk {chunk_idx}: ID={block_id}, bytes={start_byte}-{end_byte} ({chunk_mb:.2f} MB), isLast={is_last}", file=sys.stderr)

        try:
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={start_byte}-{end_byte}'
            )

            # Read entire chunk into memory
            chunk_data = response['Body'].read()
            chunk_size = len(chunk_data)

            if debug:
                print(f"[STREAM] Downloaded {chunk_size} bytes for block {block_id}", file=sys.stderr)

            # Write header (24 bytes) using existing function
            write_block_header(fifo_handle, block_id, chunk_size, is_last)
            total_bytes_written += 24

            # Write chunk data
            fifo_handle.write(chunk_data)
            fifo_handle.flush()
            total_bytes_written += chunk_size

            if debug:
                print(f"[STREAM] Wrote block {block_id}: 24-byte header + {chunk_size} data bytes", file=sys.stderr)

        except Exception as e:
            print(f"[STREAM] ERROR downloading chunk {chunk_idx} (block_id={block_id}): {e}", file=sys.stderr)
            raise

    if debug:
        print(f"[STREAM] Complete: {total_bytes_written} bytes total ({len(chunks)} chunks)", file=sys.stderr)

    return total_bytes_written


def upload_to_s3(bucket, key, data, debug=False):
    """
    Upload data to S3.
    Handles errors gracefully - logs but doesn't crash.
    """
    try:
        if debug:
            print(f"[S3_UPLOAD] Uploading to s3://{bucket}/{key} ({len(data)} bytes)", file=sys.stderr)

        s3 = boto3.client('s3')
        s3.put_object(Bucket=bucket, Key=key, Body=data)

        if debug:
            print(f"[S3_UPLOAD] Upload complete", file=sys.stderr)

        return True
    except Exception as e:
        print(f"[S3_UPLOAD] Error uploading to S3: {e}", file=sys.stderr)
        return False

import struct

def write_block_header(fifo_handle, block_id, block_size, is_last):
    """
    Write binary block header compatible with runtime/r_merge.c

    Format matches C struct with padding:
      typedef struct block_header {
        int64_t id;        // offset 0, 8 bytes
        size_t blockSize;  // offset 8, 8 bytes
        int8_t isLast;     // offset 16, 1 byte
        // 7 bytes padding   // offset 17-23
      } block_header;

    Total: 24 bytes (C struct alignment)
    """
    # Pack with 7 padding bytes: '<qQb7x'
    # '<' = little-endian
    # 'q' = int64_t (8 bytes signed)
    # 'Q' = size_t (8 bytes unsigned)
    # 'b' = int8_t (1 byte signed)
    # '7x' = 7 padding bytes
    header = struct.pack('<qQb7x', block_id, block_size, is_last)
    fifo_handle.write(header)




def main():
    try:
        # Log initialization
        log_timing("INIT", "Lambda initialization", debug=True)

        # Parse command line arguments
        if len(sys.argv) < 4:
            print("Usage: python3.9 s3-shard-reader.py <s3_key> <output_fifo> <byte_range> split=<N> num_shards=<N> job_uid=<UID>", file=sys.stderr)
            print("Example: python3.9 s3-shard-reader.py 'oneliners/inputs/1M.txt' /tmp/fifo21 bytes=0-524287 split=0 num_shards=4 job_uid=test-123", file=sys.stderr)
            sys.exit(1)

        # Positional arguments
        s3_key = sys.argv[1]
        output_fifo = sys.argv[2]
        byte_range_str = sys.argv[3]

        # Parse keyword arguments
        kwargs = parse_keyword_args(sys.argv[4:])
        shard_index = int(kwargs.get('shard', 0))
        num_shards = int(kwargs.get('num_shards', 1))
        job_uid = kwargs.get('job_uid', 'default')

        # PERF_MODE disables debug output for clean performance tests
        PERF_MODE = os.environ.get('PASH_PERF_MODE', 'false').lower() == 'true'
        debug = kwargs.get('debug', 'false').lower() == 'true' and not PERF_MODE

        # Get bucket from environment
        s3_bucket = os.environ.get('AWS_BUCKET')
        if not s3_bucket:
            print("Error: AWS_BUCKET environment variable not set", file=sys.stderr)
            sys.stderr.flush()
            sys.exit(1)

        # Parse byte ranges: single string OR JSON array
        chunks = []
        try:
            chunks_data = json.loads(byte_range_str)

            # Multi-chunk mode: JSON array
            if isinstance(chunks_data, list) and len(chunks_data) > 0:
                chunks = chunks_data

                if debug:
                    chunk_ids = [c['block_id'] for c in chunks]
                    total_mb = sum((c['end'] - c['start'] + 1) / (1024*1024) for c in chunks)
                    print(f"[MAIN {shard_index}] Multi-chunk mode: {len(chunks)} chunks, IDs {chunk_ids}, {total_mb:.2f} MB", file=sys.stderr)
            else:
                raise ValueError("Not a valid chunk list")

        except (json.JSONDecodeError, ValueError, KeyError, TypeError):
            # Single-chunk mode: backward compatible
            start_byte, end_byte = parse_byte_range(byte_range_str)
            chunks = [{
                "start": start_byte,
                "end": end_byte,
                "block_id": 0,
                "skip_first_line": False
            }]

            if debug:
                print(f"[MAIN {shard_index}] Single-chunk mode: bytes={start_byte}-{end_byte}", file=sys.stderr)

        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Starting", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] split: {shard_index}/{num_shards}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] UID: {job_uid}", file=sys.stderr)
            sys.stderr.flush()

        # Smart boundary mode: boundaries are line-aligned, no inter-lambda communication needed
        if debug and shard_index > 0:
            print(f"[{_now_ts()}][MAIN {shard_index}] Smart boundary mode: Boundary is line-aligned, skipping request listener", file=sys.stderr)
            sys.stderr.flush()

        log_timing("ARGS_PARSED", f"Parsed {len(chunks)} chunks", debug)

        # 3. Open FIFO and stream S3 directly to it
        log_timing("FIFO_OPEN_START", f"Opening FIFO {output_fifo}", debug)
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Opening FIFO: {output_fifo}", file=sys.stderr)
            sys.stderr.flush()

        # Open FIFO for writing (blocks until downstream reader connects)
        with open(output_fifo, 'wb', buffering=0) as fifo:
            log_timing("FIFO_OPEN_END", "FIFO connected to downstream", debug)
            if debug:
                print(f"[{_now_ts()}][MAIN {shard_index}] FIFO connected, streaming chunks", file=sys.stderr)
                sys.stderr.flush()

            # Stream S3 data to FIFO
            log_timing("STREAM_START", "Starting S3→FIFO stream", debug)

            bytes_streamed = stream_s3_to_fifo_handle(
                fifo, s3_bucket, s3_key, chunks,  # Pass chunks list instead of start/end
                debug=debug
            )
            log_timing("STREAM_END", f"Streamed {bytes_streamed} bytes", debug)

            if debug:
                print(f"[{_now_ts()}][MAIN {shard_index}] Streamed {bytes_streamed} bytes from S3", file=sys.stderr)
                sys.stderr.flush()

            # Smart boundary mode: boundary already line-aligned, no tail bytes needed
            if debug and shard_index < num_shards - 1 and bytes_streamed > 0:
                print(f"[{_now_ts()}][MAIN {shard_index}] Smart boundary mode: Boundary is line-aligned, skipping tail byte request", file=sys.stderr)
                sys.stderr.flush()

        # Note: S3 upload skipped in streaming mode (data already consumed by downstream)
        log_timing("FIFO_CLOSE", "Closing FIFO", debug)
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] FIFO closed, data transmitted to downstream", file=sys.stderr)
            sys.stderr.flush()

        log_timing("COMPLETE", "Lambda complete", debug)
        print_timing_summary(debug)
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Complete", file=sys.stderr)
            sys.stderr.flush()

    except Exception as e:
        print(f"[{_now_ts()}][ERROR] Fatal exception in main: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
        sys.exit(1)


if __name__ == '__main__':
    main()