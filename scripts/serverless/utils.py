import os
import boto3
import sys
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
    while True:
        response = (logs_client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
        if not response["logStreams"] or len(response["logStreams"]) == 0:
            break
        for log_stream in response["logStreams"]:
            file_name = ""
            with open(f"{out_dir}/temp.log", "w") as f:
                log_stream_response = logs_client.get_log_events(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
                file_name = f"log_id_{log_stream['logStreamName'].split(']')[1]}"
                for log in log_stream_response["events"]:
                    f.write(f"[{datetime.fromtimestamp(log['timestamp']/1000.0)}] {log['message']}")
                    if "Executing script ID " in log['message']:
                        script_id = log['message'].split("Executing script ID ")[1].split()[0]
                        if script_id:
                            file_name += f"_script_id_{script_id}"
            os.rename(f"{out_dir}/temp.log", f"{out_dir}/{file_name}.log")
            logs_client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])

def save_then_delete_scripts(out_dir: str = "scripts"):
    """
    Saves all scripts in the specified directory to files.
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    paginator = s3_client.get_paginator('list_objects_v2')
    page_iterator = paginator.paginate(Bucket='yizheng', Prefix='sls-scripts/')
    for page in page_iterator:
        if 'Contents' not in page:
            continue
        for obj in page['Contents']:
            key = obj['Key']
            print(f"Downloading {key}...")
            file_name = os.path.join(out_dir, key.split('/')[-2], key.split('/')[-1])
            file_dir = os.path.dirname(file_name)
            if not os.path.exists(file_dir):
                os.makedirs(file_dir)
            if not os.path.exists(file_name):
                s3_client.download_file('yizheng', key, file_name)
            s3_client.delete_object(Bucket='yizheng', Key=key)
            print(f"Deleted {key} from S3.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        debug_dir_prefix = sys.argv[1]
    else:
        debug_dir_prefix = "debug"
    if not os.path.exists(debug_dir_prefix):
        os.mkdir(debug_dir_prefix)
    else:
        print(f"Debug directory already exists: `{debug_dir_prefix}`, skipping everything to avoid overwriting.")
        sys.exit(0)

    logs_folder = os.path.join(debug_dir_prefix, "logs")
    scripts_folder = os.path.join(debug_dir_prefix, "scripts")
    save_then_delete_log_streams(logs_folder)
    save_then_delete_scripts(scripts_folder)