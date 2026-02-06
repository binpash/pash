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
from io import BytesIO


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


def download_s3_range(bucket, key, start_byte, end_byte):
    """
    Download a specific byte range from S3.
    Returns bytes.
    """
    s3 = boto3.client('s3')

    response = s3.get_object(
        Bucket=bucket,
        Key=key,
        Range=f'bytes={start_byte}-{end_byte}'
    )

    return response['Body'].read()


def skip_until_newline(data):
    """
    Skip bytes until first newline (inclusive).
    Port from Go client/dfs_split_reader.go lines 96-101.

    Returns the data after the first newline, or empty bytes if no newline found.
    """
    try:
        newline_pos = data.index(b'\n')
        # Return data after the newline (skip the newline itself)
        return data[newline_pos + 1:]
    except ValueError:
        # No newline found, entire chunk is skipped
        return b''


def stream_s3_bytes(bucket, key, start_byte, chunk_size=4096):
    """
    Generator that yields bytes from S3 starting at start_byte.
    Reads in chunks for efficiency but yields byte-by-byte.
    """
    s3 = boto3.client('s3')
    current_byte = start_byte

    while True:
        try:
            end_byte = current_byte + chunk_size - 1
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={current_byte}-{end_byte}'
            )
            chunk = response['Body'].read()

            if not chunk:
                break

            for byte in chunk:
                yield byte

            current_byte = end_byte + 1

        except Exception as e:
            # Hit end of file or error
            if 'InvalidRange' in str(e):
                break
            raise


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


def start_server_thread(shard_index, job_uid, s3_bucket, s3_key, byte_range_str, debug=False):
    """
    Background thread that serves process i-1's request for tail bytes.
    Uses pashlib send to transmit bytes until newline.

    Port from Go server logic (server/server.go).
    """
    def server():
        try:
            # UID for communication from i-1 to i
            comm_uid = f"{job_uid}-split-{shard_index-1}-to-{shard_index}"

            # Create send FIFO
            send_fifo = f"/tmp/pash_send_{shard_index}_{os.getpid()}"

            # Clean up any existing FIFO
            try:
                os.unlink(send_fifo)
            except FileNotFoundError:
                pass

            os.mkfifo(send_fifo)

            if debug:
                print(f"[SERVER {shard_index}] Created FIFO: {send_fifo}", file=sys.stderr)
                print(f"[SERVER {shard_index}] UID: {comm_uid}", file=sys.stderr)

            # Start pashlib send (blocks until recv side connects)
            pashlib_cmd = [PASHLIB_PATH, f"send*{comm_uid}*0*1*{send_fifo}"]
            if debug:
                print(f"[SERVER {shard_index}] Starting: {' '.join(pashlib_cmd)}", file=sys.stderr)

            proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)

            # Stream bytes from S3 to send_fifo until newline
            start_byte = parse_byte_range(byte_range_str)[0]

            if debug:
                print(f"[SERVER {shard_index}] Streaming from S3 byte {start_byte}...", file=sys.stderr)

            with open(send_fifo, 'wb') as f:
                for byte in stream_s3_bytes(s3_bucket, s3_key, start_byte):
                    f.write(bytes([byte]))
                    if byte == ord('\n'):
                        if debug:
                            print(f"[SERVER {shard_index}] Found newline, stopping", file=sys.stderr)
                        break

            # Cleanup
            proc.wait()
            os.unlink(send_fifo)

            if debug:
                print(f"[SERVER {shard_index}] Complete", file=sys.stderr)

        except Exception as e:
            print(f"[SERVER {shard_index}] Error: {e}", file=sys.stderr)
            raise

    # Run server in background thread
    thread = threading.Thread(target=server, daemon=True)
    thread.start()

    return thread


def request_tail_bytes_pashlib(shard_index, job_uid, debug=False):
    """
    Request bytes from next process until newline using pashlib recv.
    Port from Go lines 62-82 (byte-by-byte reading).
    """
    comm_uid = f"{job_uid}-split-{shard_index}-to-{shard_index+1}"
    # we want this to read from (i+1), i+2, ... n lambdas until newline

    # Create recv FIFO
    recv_fifo = f"/tmp/pash_recv_{shard_index}_{os.getpid()}"

    # Clean up any existing FIFO
    try:
        os.unlink(recv_fifo)
    except FileNotFoundError:
        pass

    os.mkfifo(recv_fifo)

    if debug:
        print(f"[CLIENT {shard_index}] Created FIFO: {recv_fifo}", file=sys.stderr)
        print(f"[CLIENT {shard_index}] UID: {comm_uid}", file=sys.stderr)

    try:
        # Start pashlib recv
        pashlib_cmd = [PASHLIB_PATH, f"recv*{comm_uid}*1*0*{recv_fifo}"]
        if debug:
            print(f"[CLIENT {shard_index}] Starting: {' '.join(pashlib_cmd)}", file=sys.stderr)

        proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)

        # Read from recv_fifo byte-by-byte until newline
        tail_bytes = bytearray()

        if debug:
            print(f"[CLIENT {shard_index}] Reading tail bytes...", file=sys.stderr)

        with open(recv_fifo, 'rb') as f: #reads byte by byte from FIFO (opt/pashlib)
            while True:
                byte = f.read(1)
                if not byte:  # EOF
                    if debug:
                        print(f"[CLIENT {shard_index}] EOF reached", file=sys.stderr)
                    break
                tail_bytes.append(byte[0])
                if byte == b'\n':
                    if debug:
                        print(f"[CLIENT {shard_index}] Found newline", file=sys.stderr)
                    break

                    #TODO. how do we tell the server to stop sending more data?
        proc.wait()

        if debug:
            print(f"[CLIENT {shard_index}] Received {len(tail_bytes)} bytes", file=sys.stderr)

        return bytes(tail_bytes)

    finally:
        # Cleanup
        try:
            os.unlink(recv_fifo)
        except FileNotFoundError:
            pass


def main():
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
    debug = kwargs.get('debug', 'false').lower() == 'true'

    # Get bucket from environment
    s3_bucket = os.environ.get('AWS_BUCKET')
    if not s3_bucket:
        print("Error: AWS_BUCKET environment variable not set", file=sys.stderr)
        sys.exit(1)

    if debug:
        print(f"[MAIN {shard_index}] Starting", file=sys.stderr)
        print(f"[MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
        print(f"[MAIN {shard_index}] Range: {byte_range_str}", file=sys.stderr)
        print(f"[MAIN {shard_index}] split: {shard_index}/{num_shards}", file=sys.stderr)
        print(f"[MAIN {shard_index}] UID: {job_uid}", file=sys.stderr)

    # 1. Start server thread (if split > 0)
    server_thread = None
    if shard_index > 0:
        if debug:
            print(f"[MAIN {shard_index}] Starting server thread", file=sys.stderr)
        server_thread = start_server_thread(
            shard_index, job_uid, s3_bucket, s3_key, byte_range_str, debug
        )

    # 2. Download S3 chunk
    start_byte, end_byte = parse_byte_range(byte_range_str)
    if debug:
        print(f"[MAIN {shard_index}] Downloading bytes {start_byte}-{end_byte}", file=sys.stderr)

    my_data = download_s3_range(s3_bucket, s3_key, start_byte, end_byte)

    if debug:
        print(f"[MAIN {shard_index}] Downloaded {len(my_data)} bytes", file=sys.stderr)

    # 3. Skip first line (if split > 0)
    if shard_index > 0:
        original_len = len(my_data)
        my_data = skip_until_newline(my_data)
        if debug:
            print(f"[MAIN {shard_index}] Skipped first line: {original_len - len(my_data)} bytes", file=sys.stderr)

    # 4. Build output buffer
    output_buffer = BytesIO()
    output_buffer.write(my_data)

    # 5. Request tail bytes via pashlib (if needed)
    if shard_index < num_shards - 1 and not my_data.endswith(b'\n'):
        if debug:
            print(f"[MAIN {shard_index}] Requesting tail bytes from split {shard_index+1}", file=sys.stderr)

        tail_bytes = request_tail_bytes_pashlib(shard_index, job_uid, debug)
        output_buffer.write(tail_bytes)

        if debug:
            print(f"[MAIN {shard_index}] Appended {len(tail_bytes)} tail bytes", file=sys.stderr)

    # 6. Write to output FIFO (BLOCKS until downstream consumer reads)
    if debug:
        print(f"[MAIN {shard_index}] Writing {output_buffer.tell()} bytes to {output_fifo}", file=sys.stderr)

    output_data = output_buffer.getvalue()

    # Normalize line endings once and for all
    # Ensures deterministic sort order across shards
    output_data = output_data.replace(b'\r\n', b'\n')

    with open(output_fifo, 'wb') as f:
        f.write(output_data)

    # 7. Upload output to S3
    # Extract base filename (e.g., "oneliners/inputs/1M.txt" -> "1M.txt")
    base_filename = os.path.basename(s3_key)
    # Split into stem and extension (e.g., "1M.txt" -> "1M", ".txt")
    stem, ext = os.path.splitext(base_filename)
    # Create output key: outputs/1M-shard-1.txt
    output_s3_key = f"outputs/{stem}-shard-{shard_index}{ext}"

    if debug:
        print(f"[MAIN {shard_index}] Uploading to S3: {output_s3_key}", file=sys.stderr)

    upload_to_s3(s3_bucket, output_s3_key, output_data, debug)

    # Wait for server thread to complete (if started)
    if server_thread:
        server_thread.join(timeout=5)

    if debug:
        print(f"[MAIN {shard_index}] Complete", file=sys.stderr)


if __name__ == '__main__':
    main()
