import os
import boto3
import sys
import collections
import shutil
import re
from datetime import datetime

logs_client = boto3.client('logs')
s3_client = boto3.client('s3')

def debug_print(enabled: bool, message: str):
    if enabled:
        print(message)

def parse_cli_args(argv):
    """
    Parse CLI arguments while preserving backward-compatible positional args:
      1) output dir prefix
      2) start_time_ms
      3) job_id
    Additional optional flags:
      --utils-logs / --util-logs / --debug-logs: enable verbose debug logs
    """
    debug_logs = False
    positional_args = []
    debug_flags = {"--utils-logs", "--util-logs", "--debug-logs"}

    for arg in argv:
        if arg in debug_flags:
            debug_logs = True
        else:
            positional_args.append(arg)

    return positional_args, debug_logs

def delete_log_streams():
    """
    Deletes all log streams in the specified log group.
    """
    while True:
        response = (logs_client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
        if not response["logStreams"] or len(response["logStreams"]) == 0:
            break
        for log_stream in response["logStreams"]:
            logs_client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])

def save_then_delete_log_streams(out_dir: str = "logs", start_time_ms: int = None, job_id: str = None, debug_logs: bool = False):
    """
    Saves all log streams in the specified log group to files.
    """
    log_group_name = "/aws/lambda/lambda"

    # Add verification output
    print(f"[Analysis] Log group: {log_group_name}")
    print(f"[Analysis] Output directory: {out_dir}")
    print(f"[Analysis] Start time filter: {start_time_ms} ({datetime.fromtimestamp(start_time_ms/1000.0) if start_time_ms else 'No filter'})")
    print(f"[Analysis] Job ID filter: {job_id if job_id else 'No filter (save all logs)'}")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    log_stream_count = 0
    total_streams_found = 0
    # streams_skipped_by_timestamp = 0
    streams_skipped_by_job_id = 0
    streams_with_no_events = 0
    
    # Outer loop: List Log Streams (Paginated)
    next_token_describe = None
    while True:
        describe_args = {
            'logGroupName': log_group_name,
            'orderBy': 'LastEventTime',
            'descending': True
        }
        if next_token_describe:
            describe_args['nextToken'] = next_token_describe

        try:
            response = logs_client.describe_log_streams(**describe_args)
        except logs_client.exceptions.ResourceNotFoundException:
            print(f"[Error] Log group {log_group_name} not found.")
            return

        if not response.get("logStreams"):
            debug_print(debug_logs, "[DEBUG] No more log streams found (pagination complete)")
            break

        page_streams = response.get("logStreams", [])
        total_streams_found += len(page_streams)
        debug_print(debug_logs, f"[DEBUG] Found {len(page_streams)} log streams in this page")

        for log_stream in page_streams:
            stream_name = log_stream["logStreamName"]

            # Skip streams from before start_time_ms (performance optimization)
            # This prevents checking hundreds of old log streams from previous runs
            # if start_time_ms and log_stream.get('lastEventTimestamp', 0) < start_time_ms:
                # streams_skipped_by_timestamp += 1
                # continue

            debug_print(debug_logs, f"[DEBUG] Checking stream: {stream_name}")
            
            # INNER LOOP: Fetch Log Events (Paginated) -> BUG FIX 1
            # NO TIMESTAMP FILTERING - fetch all events from this stream
            all_events = []
            next_token_events = None

            while True:
                get_events_args = {
                    'logGroupName': log_group_name,
                    'logStreamName': stream_name,
                    'startFromHead': True # Read from start to get the Job ID
                }
                # NO startTime filter - removed to rely solely on job_id filtering
                if next_token_events:
                    get_events_args['nextToken'] = next_token_events

                events_response = logs_client.get_log_events(**get_events_args)
                events_batch = events_response.get("events", [])

                if not events_batch:
                    break

                all_events.extend(events_batch)

                # Check if we reached the end of the stream
                if events_response['nextForwardToken'] == next_token_events:
                    break
                next_token_events = events_response['nextForwardToken']

            debug_print(debug_logs, f"[DEBUG]   Total events fetched: {len(all_events)}")

            if not all_events:
                debug_print(debug_logs, "[DEBUG]   SKIPPED: No events in stream")
                streams_with_no_events += 1
                continue

            # Check for Job ID match
            file_name = f"log_id_{stream_name.split(']')[1]}"
            file_suffix = ""
            has_matching_logs = False
            job_id_pattern = re.compile(rf'\[JOB:{re.escape(job_id)}\]') if job_id else None

            # First pass: Check match and find script ID
            for log in all_events:
                msg = log['message']
                if job_id_pattern and job_id_pattern.search(msg):
                    has_matching_logs = True
                
                if "Executing script ID " in msg:
                    try:
                        script_id = msg.split("Executing script ID ")[1].split()[0]
                        if script_id and file_suffix == "":
                            file_suffix = f"_{script_id}"
                    except IndexError:
                        pass

            # Save logs only if they match the Job ID (or if no Job ID was provided)
            if not job_id or has_matching_logs:
                debug_print(debug_logs, "[DEBUG]   MATCH: Saving log stream")
                final_path = f"{out_dir}/{file_name}{file_suffix}.log"
                with open(final_path, "w") as f:
                    for log in all_events:
                        timestamp = datetime.fromtimestamp(log['timestamp']/1000.0)
                        # BUG FIX 2: Added newline character \n
                        f.write(f"[{timestamp}] {log['message'].rstrip()}\n")
                log_stream_count += 1
            else:
                debug_print(debug_logs, f"[DEBUG]   SKIPPED SAVE: Job ID pattern '{job_id}' not found in logs")
                streams_skipped_by_job_id += 1

            # Always delete the log stream from CloudWatch once it has been processed
            logs_client.delete_log_stream(
                logGroupName=log_group_name,
                logStreamName=stream_name
            )

        # Handle pagination for describe_log_streams
        next_token_describe = response.get('nextToken')
        if not next_token_describe:
            break

    if debug_logs:
        print(f"\n[DEBUG SUMMARY]")
        print(f"  Total streams found: {total_streams_found}")
        # print(f"  Streams checked (after timestamp filter): {total_streams_found - streams_skipped_by_timestamp}")
        # print(f"  Skipped by timestamp (before {datetime.fromtimestamp(start_time_ms/1000.0) if start_time_ms else 'N/A'}): {streams_skipped_by_timestamp}")
        print(f"  Streams saved: {log_stream_count}")
        print(f"  Skipped by job ID filter: {streams_skipped_by_job_id}")
        print(f"  Skipped due to no events: {streams_with_no_events}")
    print(f"[Analysis] Total log streams saved: {log_stream_count}")
    

def save_then_delete_log_streams_copy(out_dir: str = "logs", start_time_ms: int = None, job_id: str = None):
    """
    Saves all log streams in the specified log group to files.
    Args:
        out_dir: Directory to save logs
        start_time_ms: Only fetch logs after this timestamp (milliseconds since epoch)
        job_id: Only fetch logs matching this job ID (logs with [JOB:{job_id}] prefix)
    """
    log_group_name = "/aws/lambda/lambda"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    log_stream_count = 0
    next_token = None

    while True:
        describe_args = {
            'logGroupName': log_group_name,
            'orderBy': 'LastEventTime',
            'descending': True
        }
        if next_token:
            describe_args['nextToken'] = next_token

        response = logs_client.describe_log_streams(**describe_args)

        if not response["logStreams"] or len(response["logStreams"]) == 0:
            break

        for log_stream in response["logStreams"]:
            # Skip logs from before start_time_ms
            # if start_time_ms and log_stream.get('lastEventTimestamp', 0) < start_time_ms:
                # continue

            # Add startTime parameter to filter events
            get_events_args = {
                'logGroupName': log_group_name,
                'logStreamName': log_stream["logStreamName"]
            }
            if start_time_ms:
                get_events_args['startTime'] = start_time_ms

            log_stream_response = logs_client.get_log_events(**get_events_args)

            file_name = f"log_id_{log_stream['logStreamName'].split(']')[1]}"
            file_suffix = ""

            # If job_id is specified, first check if this log stream contains matching logs
            has_matching_logs = False
            job_id_pattern = re.compile(rf'\[JOB:{re.escape(job_id)}\]') if job_id else None

            # First pass: determine if this log stream belongs to this job
            if job_id_pattern:
                for log in log_stream_response["events"]:
                    if job_id_pattern.search(log['message']):
                        has_matching_logs = True
                        break

            # Second pass: if no job_id filter OR if this stream matches, write all logs
            if not job_id or has_matching_logs:
                with open(f"{out_dir}/temp.log", "w") as f:
                    for log in log_stream_response["events"]:
                        f.write(f"[{datetime.fromtimestamp(log['timestamp']/1000.0)}] {log['message']}")

                        # Extract script ID for filename
                        if "Executing script ID " in log['message']:
                            script_id = log['message'].split("Executing script ID ")[1].split()[0]
                            if script_id and file_suffix == "":
                                file_suffix = f"_{script_id}"

                # Save and delete the log stream
                os.rename(f"{out_dir}/temp.log", f"{out_dir}/{file_name}{file_suffix}.log")
                logs_client.delete_log_stream(
                    logGroupName=log_group_name,
                    logStreamName=log_stream["logStreamName"]
                )
                log_stream_count += 1
            else: #stream doesn't match job_id, skip it and don't delete it from CloudWatch
                logs_client.delete_log_stream(
                    logGroupName=log_group_name,
                    logStreamName=log_stream["logStreamName"]
                )

        # Check if there are more pages
        next_token = response.get('nextToken')
        if not next_token:
            break  # No more pages, exit loop

    print(f"[Analysis] Total log streams saved: {log_stream_count}")

def save_then_delete_scripts(out_dir: str = "scripts"):
    """
    Saves all scripts in the specified directory to files.
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    script_count = 0
    paginator = s3_client.get_paginator('list_objects_v2')
    page_iterator = paginator.paginate(Bucket=os.environ['AWS_BUCKET'], Prefix='sls-scripts/')
    for page in page_iterator:
        if 'Contents' not in page:
            continue
        for obj in page['Contents']:
            key = obj['Key']
            file_name = os.path.join(out_dir, key.split('/')[-2], key.split('/')[-1])
            file_dir = os.path.dirname(file_name)
            if not os.path.exists(file_dir):
                os.makedirs(file_dir)
            if not os.path.exists(file_name):
                s3_client.download_file(os.environ['AWS_BUCKET'], key, file_name)
                script_count += 1
            s3_client.delete_object(Bucket=os.environ['AWS_BUCKET'], Key=key)
    print(f"[Analysis] Total shell scripts saved: {script_count}")

def analyze_logs(logs_folder: str):
    """
    Analyzes the logs in the specified folder and prints the summary.
    """
    if not os.path.exists(logs_folder):
        print(f"[Analysis] Logs folder does not exist: {logs_folder}")
        return

    total_billed_time = []
    count_billed_time = 0
    err_log_files = []
    err_type = collections.defaultdict(str)

    print(f"[Analysis] Analyzing logs in folder: {logs_folder}")

    for log_file in os.listdir(logs_folder):
        err_found = False
        with open(os.path.join(logs_folder, log_file), 'r') as f:
            for line in f:
                # Broad detection for errors and process termination signals
                if ("Err" in line or "err" in line or "panic" in line or
                    "Killed" in line or "Segmentation" in line or "Terminated" in line or
                    "Aborted" in line or "Abort" in line or
                    "Floating point exception" in line or
                    "Illegal instruction" in line or "Bus error" in line or
                    "Signal" in line or "Hangup" in line or
                    "core dumped" in line):
                    err_found = True

                # Set specific error types for categorization
                if "Connection reset by peer" in line:
                    err_type[log_file] = "Connection reset by peer"
                if "already in use" in line or "Address already in use" in line:
                    err_type[log_file] = "Port already in use"
                if "Timeout" in line or "timed out" in line.lower():
                    err_type[log_file] = "Timeout"
                if "OutOfMemory" in line or "out of memory" in line.lower():
                    err_type[log_file] = "Out of memory"
                if "thread 'tokio-runtime-worker' panicked" in line or "thread 'main' panicked" in line:
                    err_type[log_file] = "Thread panic"
                if line.strip().startswith("ERROR") or ": error:" in line.lower():
                    err_type[log_file] = "Lambda error"
                if "Task timed out" in line:
                    err_type[log_file] = "Task timeout"
                if "JoinError" in line:
                    err_type[log_file] = "Join error"

                # Process termination signals
                if "Killed" in line:
                    err_type[log_file] = "Process killed (OOM/Signal)"
                if "Segmentation fault" in line or "Segmentation" in line:
                    err_type[log_file] = "Segmentation fault"
                if "Terminated" in line and "Timeout" not in line:
                    err_type[log_file] = "Process terminated (Signal)"
                if "Aborted" in line or "Abort" in line:
                    err_type[log_file] = "Process aborted"
                if "Floating point exception" in line:
                    err_type[log_file] = "Floating point exception"
                if "Illegal instruction" in line:
                    err_type[log_file] = "Illegal instruction"
                if "Bus error" in line:
                    err_type[log_file] = "Bus error"
                if "core dumped" in line:
                    err_type[log_file] = "Core dumped"

                # Keep billed duration tracking
                if "Billed Duration:" in line:
                    try:
                        billed_time = int(line.split("Billed Duration:")[1].split(" ms")[0])
                        total_billed_time.append(billed_time)
                        count_billed_time += 1
                        # print(f"Billed time: {log_file}: {total_billed_time[-1]} ms")
                    except ValueError:
                        print(f"[Analysis] Could not parse billed time from line: {line.strip()}")

        if err_found:
            err_log_files.append(log_file)
    print(f"[Analysis] Total billed time: {sum(total_billed_time)} ms")
    print(f"[Analysis] Total number of billed time entries: {count_billed_time}")
    print(f"[Analysis] Cost estimate: ${sum(total_billed_time) * 0.000000028849902:.6f}")
    print(f"[Analysis] Errors may be found in the following log files:")
    for err_log_file in err_log_files:
        print(f"  {err_type[err_log_file]}:")
        print(f"    {logs_folder}/{err_log_file}")

if __name__ == "__main__":
    positional_args, debug_logs = parse_cli_args(sys.argv[1:])

    if len(positional_args) > 0:
        debug_dir_prefix = positional_args[0]
    else:
        debug_dir_prefix = "debug"

    # Read start time from command line
    start_time_ms = None
    if len(positional_args) > 1:
        try:
            start_time_ms = int(positional_args[1])
            print(f"[Analysis] Filtering logs from timestamp: {start_time_ms} ({datetime.fromtimestamp(start_time_ms/1000.0)})")
        except ValueError:
            print(f"[Analysis] Warning: Could not parse start_time_ms from '{positional_args[1]}', fetching all logs")
            start_time_ms = None

    # Read job_id from command line
    job_id = None
    if len(positional_args) > 2:
        job_id = positional_args[2]
        print(f"[Analysis] Filtering logs by job ID: {job_id}")

    # Create base directory
    os.makedirs(debug_dir_prefix, exist_ok=True)

    # Clear subdirectories to ensure fresh data
    logs_folder = os.path.join(debug_dir_prefix, "logs")
    scripts_folder = os.path.join(debug_dir_prefix, "scripts")

    if os.path.exists(logs_folder):
        shutil.rmtree(logs_folder)
    if os.path.exists(scripts_folder):
        shutil.rmtree(scripts_folder)

    # Now fetch fresh data with time and job_id filtering
    save_then_delete_log_streams(logs_folder, start_time_ms, job_id, debug_logs=debug_logs)
    save_then_delete_scripts(scripts_folder)
    analyze_logs(logs_folder)
