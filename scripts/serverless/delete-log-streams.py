import os
import boto3
from datetime import datetime
client = boto3.client('logs')
response = (client.describe_log_streams(logGroupName="/aws/lambda/lambda"))
for log_stream in response["logStreams"]:
    client.delete_log_stream(logGroupName="/aws/lambda/lambda", logStreamName=log_stream["logStreamName"])
