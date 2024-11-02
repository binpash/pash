import boto3
import sys
import json
import os
import time

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
BUCKET=os.environ.get("AWS_BUCKET")
QUEUE=os.environ.get("AWS_QUEUE")
object_key, infile, s3_folder_id = sys.argv[1:]

print(f"[s3-put-object.py] Start putting {object_key}")
start_time = time.time()

session = boto3.Session()
s3 = session.client("s3")


# Prepare the data to be stored in the S3 object
with open(infile, "rb") as file:
    object_data = file.read()

# Upload the object to S3
s3.put_object(Bucket=BUCKET, Key=object_key, Body=object_data)

end_time = time.time()
total_time = end_time - start_time
total_time = round(total_time, 3)

print(f"[s3-put-object.py] Finish putting {object_key} in {total_time} seconds")

# notify the main shell that job is done
# sqs_client = boto3.client("sqs")
# message_body = {"message": "done","output_file_id":object_key}
# try:
#     response = sqs_client.send_message(
#         QueueUrl=f"https://sqs.us-east-1.amazonaws.com/{AWS_ACCOUNT_ID}/{QUEUE}", MessageBody=json.dumps(message_body)
#     )
# except Exception as e:
#     print(e)
#   

dynamodb = session.client("dynamodb")
try:
    response = dynamodb.put_item(
        TableName="sls-result",
        Item={
            "key": {"S": s3_folder_id},
            "output_id": {"S": object_key},
        },
    )
    print(f"[s3-put-object.py] Successfully wrote item {s3_folder_id}:{object_key} to dynamoDB.")
except Exception as e:
    print(f"[s3-put-object.py] Failed to write item {s3_folder_id}:{object_key} to dynamoDB. {e}")