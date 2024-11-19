import boto3
import sys
import json
import os
import time

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
BUCKET=os.environ.get("AWS_BUCKET")
QUEUE=os.environ.get("AWS_QUEUE")
object_key, infile, s3_folder_id = sys.argv[1:]

if object_key == "/dev/null":
    print(f"[s3-put-object.py] Get /dev/null as key, ignore.")
    exit(0)

s3 = boto3.client('s3')
start_time = time.time()
with open(infile, "rb") as f:
    s3.upload_fileobj(f, BUCKET, object_key)

end_time = time.time()
total_time = end_time - start_time
total_time = round(total_time, 3)

print(f"[s3-put-object.py] Finish putting {object_key} in {total_time} seconds")

# session = boto3.Session()
# dynamodb = session.client("dynamodb")
# try:
#     response = dynamodb.put_item(
#         TableName="sls-result",
#         Item={
#             "key": {"S": s3_folder_id},
#             "output_id": {"S": object_key},
#         },
#     )
#     print(f"[s3-put-object.py] Successfully wrote item {s3_folder_id}:{object_key} to dynamoDB.")
# except Exception as e:
#     print(f"[s3-put-object.py] Failed to write item {s3_folder_id}:{object_key} to dynamoDB. {e}")
