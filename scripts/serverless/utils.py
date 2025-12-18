import os
import boto3
import sys
import collections
import shutil
from datetime import datetime

logs_client = boto3.client('logs')
s3_client = boto3.client('s3')

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

def save_then_delete_log_streams(out_dir: str = "logs"):
    """
    Saves all log streams in the specified log group to files.
    """
    log_group_name = "/aws/lambda/lambda"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    log_stream_count = 0
    while True:
        response = (logs_client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
        if not response["logStreams"] or len(response["logStreams"]) == 0:
            break
        for log_stream in response["logStreams"]:
            log_stream_count += 1
            file_name = ""
            with open(f"{out_dir}/temp.log", "w") as f:
                log_stream_response = logs_client.get_log_events(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
                file_name = f"log_id_{log_stream['logStreamName'].split(']')[1]}"
                file_suffix = ""
                for log in log_stream_response["events"]:
                    f.write(f"[{datetime.fromtimestamp(log['timestamp']/1000.0)}] {log['message']}")
                    if "Executing script ID " in log['message']:
                        script_id = log['message'].split("Executing script ID ")[1].split()[0]
                        if script_id and file_suffix == "":
                            file_suffix = f"_{script_id}"
            os.rename(f"{out_dir}/temp.log", f"{out_dir}/{file_name}{file_suffix}.log")
            logs_client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
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
                if "Err" in line or "err" in line or "panic" in line or "space" in line:
                    err_found = True
                if "Connection reset by peer" in line:
                    err_type[log_file] = "Connection reset by peer"
                if "already in use" in line:
                    err_type[log_file] = "Port already in use"
                if "Timeout" in line:
                    err_type[log_file] = "Timeout"
                if "OutOfMemory" in line:
                    err_type[log_file] = "Out of memory"
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
    if len(sys.argv) > 1:
        debug_dir_prefix = sys.argv[1]
    else:
        debug_dir_prefix = "debug"

    # Create base directory
    os.makedirs(debug_dir_prefix, exist_ok=True)

    # Clear subdirectories to ensure fresh data
    logs_folder = os.path.join(debug_dir_prefix, "logs")
    scripts_folder = os.path.join(debug_dir_prefix, "scripts")

    if os.path.exists(logs_folder):
        shutil.rmtree(logs_folder)
    if os.path.exists(scripts_folder):
        shutil.rmtree(scripts_folder)

    # Now fetch fresh data
    save_then_delete_log_streams(logs_folder)
    save_then_delete_scripts(scripts_folder)
    analyze_logs(logs_folder)
