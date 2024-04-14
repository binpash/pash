import boto3
import sys
import json
import os

object_key, infile = sys.argv[1:]
print("[s3-put-object.py] Remote write",object_key)
session = boto3.Session()

s3 = session.client("s3")

# Prepare the data to be stored in the S3 object
with open(infile, "rb") as file:
    object_data = file.read()

# Upload the object to S3
s3.put_object(Bucket="bucket", Key=object_key, Body=object_data)


# notify the main shell that job is done
sqs_client = boto3.client("sqs")
message_body = {"message": "done","output_file_id":object_key}
try:
    response = sqs_client.send_message(
        QueueUrl=f"https://sqs.us-east-1.amazonaws.com/{os.environ.get('AWS_ACCOUNT_ID')}/queue ", MessageBody=json.dumps(message_body)
    )
except Exception as e:
    print(e)
