import boto3
import sys
import json

object_key, infile, last_subgraph = sys.argv[1:]
print("Remote write",object_key)
session = boto3.Session()

s3 = session.client("s3")

# Prepare the data to be stored in the S3 object
with open(infile, "rb") as file:
    object_data = file.read()

# Upload the object to S3
s3.put_object(Bucket="nikpag", Key=object_key, Body=object_data)

if last_subgraph=="1":
    # notify the main shell that job is done
    sqs_client = boto3.client("sqs")
    message_body = {"message": "done","output_file_id":object_key}
    try:
        response = sqs_client.send_message(
            QueueUrl="https://sqs.us-east-1.amazonaws.com/347768412644/queue ", MessageBody=json.dumps(message_body)
        )
    except Exception as e:
        print(e)
