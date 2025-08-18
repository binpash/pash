import boto3
import sys
import json
import os
import time

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
BUCKET=os.environ.get("AWS_BUCKET")
object_key, infile = sys.argv[1:][:2]
DEBUG=False

if object_key == "/dev/null":
    if DEBUG:
        print(f"[s3-put-object.py] Get /dev/null as key, ignore.")
    exit(0)

if DEBUG:
    print(f"[s3-put-object.py] Start putting {object_key}")
    start_time = time.time()

s3 = boto3.client('s3')
with open(infile, "rb") as f:
    s3.upload_fileobj(f, BUCKET, object_key)

if DEBUG:
    end_time = time.time()
    total_time = end_time - start_time
    total_time = round(total_time, 3)
    print(f"[s3-put-object.py] Finish putting {object_key} in {total_time} seconds")