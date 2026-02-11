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
from botocore.exceptions import ClientError
from io import BytesIO


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


def skip_until_newline(bucket, key, start_byte, end_byte, initial_data):
    """
    Skip bytes until first newline (inclusive), reading from S3 if needed.
    Port from Go client/dfs_split_reader.go lines 96-101:

        if skipFirstLine {
            _, err = reader.ReadString('\n') //discarded
            if err != nil {
                return err
            }
        }

    In the Go version, reader.ReadString('\n') reads from the file stream
    until finding '\n', potentially reading beyond the logical block boundary.

    In the S3 version with byte ranges, we must do the same: read forward from S3
    until we find the newline that ends the partial line from the previous shard.

    Args:
        bucket: S3 bucket
        key: S3 key
        start_byte: Starting byte position for this shard
        end_byte: Ending byte position for this shard
        initial_data: Already downloaded data for this shard

    Returns the data after the first newline within the shard's range.
    If the newline is found beyond our range, returns empty bytes (entire shard was part of previous line).
    """
    try:
        # First check if newline exists in already-downloaded data
        newline_pos = initial_data.index(b'\n')
        # Return data after the newline (skip the newline itself)
        return initial_data[newline_pos + 1:]
    except ValueError:
        # No newline in initial data - need to read forward from S3
        # This handles lines longer than the shard's byte range
        bytes_skipped = len(initial_data)
        current_byte = start_byte + bytes_skipped

        # Read forward until we find newline (may go beyond our shard's range!)
        # This matches Go's reader.ReadString('\n') behavior
        for byte in stream_s3_bytes(bucket, key, current_byte):
            bytes_skipped += 1

            if byte == ord('\n'):
                # Found newline! Now return remaining data within our shard's range
                remaining_start = start_byte + bytes_skipped
                if remaining_start <= end_byte:
                    # We have data remaining in our range after skipping the line
                    return download_s3_range(bucket, key, remaining_start, end_byte)
                else:
                    # Newline was beyond our range - entire shard is part of previous line
                    return b''

        # Reached EOF without finding newline - return empty
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


def start_single_server_thread(shard_index, requester_index, job_uid, s3_bucket, s3_key, byte_range_str, debug=False, request_key=None, ready_key=None):
    """
    Background thread that serves requester_index's request for tail bytes.
    Uses pashlib send to transmit bytes until newline.

    Port from Go server logic (server/server.go).
    """
    def server():
        try:
            # UID for communication from requester_index to shard_index
            comm_uid = f"{job_uid}-split-{requester_index}-to-{shard_index}"

            # Create send FIFO (unique per requester)
            send_fifo = f"/tmp/pash_send_{shard_index}_from_{requester_index}_{os.getpid()}"

            # Clean up any existing FIFO
            try:
                os.unlink(send_fifo)
            except FileNotFoundError:
                pass

            os.mkfifo(send_fifo)

            if debug:
                print(f"[{_now_ts()}][SERVER {shard_index}] Created FIFO for requester {requester_index}: {send_fifo}", file=sys.stderr)
                print(f"[{_now_ts()}][SERVER {shard_index}] UID: {comm_uid}", file=sys.stderr)

            # Start pashlib send (blocks until recv side connects)
            pashlib_cmd = [PASHLIB_PATH, f"send*{comm_uid}*0*1*{send_fifo}"]
            if debug:
                print(f"[{_now_ts()}][SERVER {shard_index}] Starting: {' '.join(pashlib_cmd)}", file=sys.stderr)

            proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)
            if ready_key:
                s3 = boto3.client('s3')
                _s3_put_marker(s3, s3_bucket, ready_key, b"ready")
                if debug:
                    print(f"[{_now_ts()}][SERVER {shard_index}] Ready marker written for requester {requester_index}", file=sys.stderr)

            # Stream bytes from S3 to send_fifo until newline
            start_byte = parse_byte_range(byte_range_str)[0]

            if debug:
                print(f"[{_now_ts()}][SERVER {shard_index}] Streaming from S3 byte {start_byte} to requester {requester_index}...", file=sys.stderr)

            fifo_fd = os.open(send_fifo, os.O_WRONLY)

            with os.fdopen(fifo_fd, 'wb', buffering=0) as f:
                for byte in stream_s3_bytes(s3_bucket, s3_key, start_byte):
                    f.write(bytes([byte]))
                    if byte == ord('\n'):
                        if debug:
                            print(f"[{_now_ts()}][SERVER {shard_index}] Found newline, stopping send to requester {requester_index}", file=sys.stderr)
                        break

            # Cleanup
            proc.wait()
            os.unlink(send_fifo)

            if debug:
                print(f"[{_now_ts()}][SERVER {shard_index}] Complete serving requester {requester_index}", file=sys.stderr)

        except Exception as e:
            print(f"[{_now_ts()}][SERVER {shard_index}] Error serving requester {requester_index}: {e}", file=sys.stderr)
            # Don't raise - let other server threads continue
        finally:
            if request_key or ready_key:
                s3 = boto3.client('s3')
                if request_key:
                    _s3_delete_marker(s3, s3_bucket, request_key)
                if ready_key:
                    _s3_delete_marker(s3, s3_bucket, ready_key)
                if debug:
                    print(f"[{_now_ts()}][SERVER {shard_index}] Cleaned rendezvous markers for requester {requester_index}", file=sys.stderr)

    # Run server in background thread (non-daemon to keep lambda alive while serving)
    thread = threading.Thread(target=server, daemon=False)
    thread.start()

    return thread


def start_request_listener(shard_index, job_uid, s3_bucket, s3_key, byte_range_str, debug=False):
    """
    Start a background listener that launches server threads only when requested.
    """
    def listener():
        s3 = boto3.client('s3')
        served = set()
        poll_secs = max(0.05, _get_env_timeout("PASH_S3_TAIL_DISCOVERY_POLL_SECS", 0.2))

        while True:
            for requester_index in range(shard_index):
                if requester_index in served:
                    continue
                req_key = _rendezvous_key(job_uid, requester_index, shard_index, "req")
                if _s3_marker_exists(s3, s3_bucket, req_key):
                    ready_key = _rendezvous_key(job_uid, requester_index, shard_index, "ready")
                    if debug:
                        print(f"[{_now_ts()}][SERVER {shard_index}] Request detected from {requester_index}", file=sys.stderr)
                    start_single_server_thread(
                        shard_index,
                        requester_index,
                        job_uid,
                        s3_bucket,
                        s3_key,
                        byte_range_str,
                        debug,
                        request_key=req_key,
                        ready_key=ready_key
                    )
                    served.add(requester_index)
            if len(served) == shard_index:
                break
            time.sleep(poll_secs)

    thread = threading.Thread(target=listener, daemon=True)
    thread.start()
    return [thread]


def request_tail_bytes_pashlib(shard_index, num_shards, job_uid, s3_bucket, debug=False):
    """
    Request bytes from subsequent processes until newline using pashlib recv.
    Strict semantics:
      - If the next shard never replies, stop immediately.
      - If it replies but has no newline, continue to the next shard.
    Port from Go lines 62-82 (byte-by-byte reading).
    """
    accumulated_bytes = bytearray()
    s3 = boto3.client('s3')
    discovery_timeout = _get_env_timeout("PASH_S3_TAIL_DISCOVERY_TIMEOUT_SECS", 10.0)
    poll_secs = max(0.05, _get_env_timeout("PASH_S3_TAIL_DISCOVERY_POLL_SECS", 0.2))

    # Loop through subsequent shards: i+1, i+2, ..., n-1
    for target_shard in range(shard_index + 1, num_shards):
        comm_uid = f"{job_uid}-split-{shard_index}-to-{target_shard}"
        req_key = _rendezvous_key(job_uid, shard_index, target_shard, "req")
        ready_key = _rendezvous_key(job_uid, shard_index, target_shard, "ready")
        _s3_put_marker(s3, s3_bucket, req_key, b"req")
        if debug:
            print(f"[{_now_ts()}][CLIENT {shard_index}] Requested tail from shard {target_shard}", file=sys.stderr)

        if not _wait_for_marker(s3, s3_bucket, ready_key, discovery_timeout, poll_secs, debug, f"[CLIENT {shard_index}]"):
            if debug:
                print(f"[{_now_ts()}][CLIENT {shard_index}] No reply from shard {target_shard}; stopping tail search", file=sys.stderr)
            _s3_delete_marker(s3, s3_bucket, req_key)
            break
        recv_fifo = f"/tmp/pash_recv_{shard_index}_to_{target_shard}_{os.getpid()}"

        # Clean up any existing FIFO
        try:
            os.unlink(recv_fifo)
        except FileNotFoundError:
            pass

        os.mkfifo(recv_fifo)

        if debug:
            print(f"[{_now_ts()}][CLIENT {shard_index}] Attempting to read from shard {target_shard}", file=sys.stderr)
            print(f"[{_now_ts()}][CLIENT {shard_index}] Created FIFO: {recv_fifo}", file=sys.stderr)
            print(f"[{_now_ts()}][CLIENT {shard_index}] UID: {comm_uid}", file=sys.stderr)

        try:
            # Start pashlib recv for this hop
            pashlib_cmd = [PASHLIB_PATH, f"recv*{comm_uid}*1*0*{recv_fifo}"]
            if debug:
                print(f"[{_now_ts()}][CLIENT {shard_index}] Starting: {' '.join(pashlib_cmd)}", file=sys.stderr)

            proc = subprocess.Popen(pashlib_cmd, stderr=subprocess.PIPE if not debug else None)

            # Read from this shard until newline OR EOF
            found_newline = False
            bytes_from_this_shard = 0

            if debug:
                print(f"[{_now_ts()}][CLIENT {shard_index}] Reading tail bytes from shard {target_shard}...", file=sys.stderr)

            recv_fd = os.open(recv_fifo, os.O_RDONLY)
            try:
                while True:
                    byte = os.read(recv_fd, 1)
                    if not byte:  # EOF from this shard
                        if debug:
                            print(f"[{_now_ts()}][CLIENT {shard_index}] EOF from shard {target_shard} after {bytes_from_this_shard} bytes", file=sys.stderr)
                        break

                    accumulated_bytes.append(byte[0])
                    bytes_from_this_shard += 1

                    if byte == b'\n':  # Found newline!
                        found_newline = True
                        if debug:
                            print(f"[{_now_ts()}][CLIENT {shard_index}] Found newline in shard {target_shard}", file=sys.stderr)
                        break
            finally:
                os.close(recv_fd)

            proc.wait()

            # If we found newline, we're done - break out of loop
            if found_newline:
                if debug:
                    print(f"[{_now_ts()}][CLIENT {shard_index}] Total received: {len(accumulated_bytes)} bytes across {target_shard - shard_index} hop(s)", file=sys.stderr)
                break

            # Otherwise, continue to next shard
            if debug:
                print(f"[{_now_ts()}][CLIENT {shard_index}] No newline in shard {target_shard}, continuing to next shard...", file=sys.stderr)

        finally:
            # Cleanup FIFO for this hop
            try:
                os.unlink(recv_fifo)
            except FileNotFoundError:
                pass
            _s3_delete_marker(s3, s3_bucket, req_key)
            _s3_delete_marker(s3, s3_bucket, ready_key)

    # After loop: either found newline OR exhausted all shards
    if not accumulated_bytes.endswith(b'\n') and debug:
        print(f"[{_now_ts()}][CLIENT {shard_index}] WARNING: No newline found in any subsequent shard", file=sys.stderr)

    return bytes(accumulated_bytes)


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
        print(f"[{_now_ts()}][MAIN {shard_index}] Starting", file=sys.stderr)
        print(f"[{_now_ts()}][MAIN {shard_index}] S3: {s3_bucket}/{s3_key}", file=sys.stderr)
        print(f"[{_now_ts()}][MAIN {shard_index}] Range: {byte_range_str}", file=sys.stderr)
        print(f"[{_now_ts()}][MAIN {shard_index}] split: {shard_index}/{num_shards}", file=sys.stderr)
        print(f"[{_now_ts()}][MAIN {shard_index}] UID: {job_uid}", file=sys.stderr)

    # 1. Start request listener (if split > 0)
    server_threads = []
    if shard_index > 0:
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Starting request listener for {shard_index} potential requesters", file=sys.stderr)
        server_threads = start_request_listener(
            shard_index, job_uid, s3_bucket, s3_key, byte_range_str, debug
        )

    # 2. Download S3 chunk
    start_byte, end_byte = parse_byte_range(byte_range_str)
    if debug:
        print(f"[{_now_ts()}][MAIN {shard_index}] Downloading bytes {start_byte}-{end_byte}", file=sys.stderr)

    my_data = download_s3_range(s3_bucket, s3_key, start_byte, end_byte)

    if debug:
        print(f"[{_now_ts()}][MAIN {shard_index}] Downloaded {len(my_data)} bytes", file=sys.stderr)

    # 3. Skip first line (if split > 0)
    if shard_index > 0:
        original_len = len(my_data)
        my_data = skip_until_newline(s3_bucket, s3_key, start_byte, end_byte, my_data)
        if debug:
            skipped_bytes = original_len - len(my_data)
            print(f"[{_now_ts()}][MAIN {shard_index}] Skipped first line: {skipped_bytes} bytes", file=sys.stderr)
            if skipped_bytes > original_len:
                print(f"[{_now_ts()}][MAIN {shard_index}] (Read {skipped_bytes - original_len} additional bytes from S3 to find newline)", file=sys.stderr)

    # 4. Build output buffer
    output_buffer = BytesIO()
    output_buffer.write(my_data)

    # 5. Request tail bytes via pashlib (if needed)
    # Only request tail bytes if we have data that doesn't end with newline
    # If my_data is empty, the entire shard was part of a previous line - don't request tail bytes
    if shard_index < num_shards - 1 and my_data:# and not my_data.endswith(b'\n'):
        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Requesting tail bytes (multi-hop from {shard_index+1} onwards)", file=sys.stderr)

        tail_bytes = request_tail_bytes_pashlib(shard_index, num_shards, job_uid, s3_bucket, debug)
        output_buffer.write(tail_bytes)

        if debug:
            print(f"[{_now_ts()}][MAIN {shard_index}] Appended {len(tail_bytes)} tail bytes", file=sys.stderr)

    # 6. Write to output FIFO (BLOCKS until downstream consumer reads)
    if debug:
        print(f"[{_now_ts()}][MAIN {shard_index}] Writing {output_buffer.tell()} bytes to {output_fifo}", file=sys.stderr)

    output_data = output_buffer.getvalue()


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
        print(f"[{_now_ts()}][MAIN {shard_index}] Uploading to S3: {output_s3_key}", file=sys.stderr)

    upload_to_s3(s3_bucket, output_s3_key, output_data, debug)

    if debug and server_threads:
        print(f"[{_now_ts()}][MAIN {shard_index}] {len(server_threads)} server threads still running (waiting for them to finish)", file=sys.stderr)

    if debug:
        print(f"[{_now_ts()}][MAIN {shard_index}] Complete", file=sys.stderr)


if __name__ == '__main__':
    main()
