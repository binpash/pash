import os
import boto3
from datetime import datetime
client = boto3.client('logs')
while True:
    response = (client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
    if not response["logStreams"] or len(response["logStreams"]) == 0:
        break
    for log_stream in response["logStreams"]:
        client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
