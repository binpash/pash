#!/usr/bin/env python3.9
"""
Phase 2: S3 split Reader with Pashlib Inter-Process Communication

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




def stream_s3_to_fifo_handle(fifo_handle, bucket, key, start_byte, end_byte, debug=False):
    """
    Stream S3 byte range DIRECTLY to the FIFO.

    FIXED: Uses s3.get_object + shutil.copyfileobj to avoid the
    'Invalid extra_args key Range' error from download_fileobj,
    while maintaining the original "pipe" behavior.

    Writes r_merge-compatible headers before each chunk.
    Only the LAST chunk gets isLast=1 to prevent data interleaving.
    """
    s3 = boto3.client('s3')

    if debug:
         print(f"[{_now_ts()}][STREAM] Strategy: DIRECT PIPE (Low Latency via get_object)", file=sys.stderr)
         sys.stderr.flush()

    # 1. Get the stream from S3 with the specific range
    # This corresponds to your original architecture's simplicity
    response = s3.get_object(
        Bucket=bucket,
        Key=key,
        Range=f'bytes={start_byte}-{end_byte}'
    )

    shutil.copyfileobj(response['Body'], fifo_handle, length=CHUNK_SIZE)
    return end_byte - start_byte + 1

    # 2. Stream data with headers using one-chunk lookahead
    # This allows us to set isLast=1 only for the final chunk
    # without buffering all data in memory
    block_id = 0
    total_bytes = 0
    prev_chunk = None

    while True:
        chunk = response['Body'].read(CHUNK_SIZE)

        if not chunk:
            break

        if prev_chunk is not None:
            # Write previous chunk with isLast=0 (not the last)
            write_block_header(fifo_handle, block_id, len(prev_chunk), 0)
            fifo_handle.write(prev_chunk)
            block_id += 1
            total_bytes += len(prev_chunk)

       
        prev_chunk = chunk

    # Write final chunk with isLast=1
    #before adding this code here we were getting lamost the right number of bytes except a bit more because it didnt know it was th elast chunk
    #now with this code were getting double 
    if prev_chunk:
        write_block_header(fifo_handle, block_id, len(prev_chunk), 1)
        fifo_handle.write(prev_chunk)
        total_bytes += len(prev_chunk)

    return total_bytes


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

        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Starting", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] Range: {byte_range_str}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] split: {shard_index}/{num_shards}", file=sys.stderr)
            print(f"[{_now_ts()}][MAIN {shard_index}] UID: {job_uid}", file=sys.stderr)
            sys.stderr.flush()

        # Smart boundary mode: boundaries are line-aligned, no inter-lambda communication needed
        if debug and shard_index > 0:
            print(f"[{_now_ts()}][MAIN {shard_index}] Smart boundary mode: Boundary is line-aligned, skipping request listener", file=sys.stderr)
            sys.stderr.flush()

        # 2. Parse byte range
        start_byte, end_byte = parse_byte_range(byte_range_str)
        log_timing("ARGS_PARSED", f"byte_range={start_byte}-{end_byte}", debug)

        # 3. Open FIFO and stream S3 directly to it
        log_timing("FIFO_OPEN_START", f"Opening FIFO {output_fifo}", debug)
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Opening FIFO: {output_fifo}", file=sys.stderr)
            sys.stderr.flush()

        # Open FIFO for writing (blocks until downstream reader connects)
        with open(output_fifo, 'wb', buffering=0) as fifo:
            log_timing("FIFO_OPEN_END", "FIFO connected to downstream", debug)
            if debug:
                print(f"[{_now_ts()}][MAIN {shard_index}] FIFO connected, streaming S3 range {start_byte}-{end_byte}", file=sys.stderr)
                sys.stderr.flush()

            # Stream S3 data to FIFO
            log_timing("STREAM_START", "Starting S3→FIFO stream", debug)
            
            bytes_streamed = stream_s3_to_fifo_handle(
                fifo, s3_bucket, s3_key, start_byte, end_byte,
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