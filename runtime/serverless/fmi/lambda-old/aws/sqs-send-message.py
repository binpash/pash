import boto3
import sys
import json
import time

queue_url, message = sys.argv[1:]

sqs_client = boto3.client("sqs")

message_body = {"message": message}

try:
    response = sqs_client.send_message(
        QueueUrl=queue_url,
        MessageBody=json.dumps(message_body),
        MessageGroupId="group",
        MessageDeduplicationId=str(time.time_ns()),
    )
except Exception as e:
    print(e)
