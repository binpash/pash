import os
import boto3
from datetime import datetime
client = boto3.client('logs')

# folder = os.getenv("LOG_FOLDER")
folder = "logs"
if not os.path.exists(folder):
    os.mkdir(folder)

while True:
    response = (client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
    if not response["logStreams"] or len(response["logStreams"]) == 0:
        break
    for log_stream in response["logStreams"]:
        file_name = ""
        with open(f"{folder}/temp.log", "w") as f:
            log_stream_response = client.get_log_events(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
            file_name = log_stream["logStreamName"].split("]")[1]
            for log in log_stream_response["events"]:
                f.write(f"[{datetime.fromtimestamp(log['timestamp']/1000.0)}] {log['message']}")
            f.write("=========================================================================\n\n\n")
        os.rename(f"{folder}/temp.log", f"{folder}/aws_lambda_log_{file_name}.log")
        client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
