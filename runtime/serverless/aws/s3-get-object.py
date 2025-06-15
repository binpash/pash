import boto3
import sys
import time
import os

BUCKET=os.environ.get("AWS_BUCKET")
object_key, outfile = sys.argv[1:]
DEBUG=False

# DEBUG
if DEBUG:
    print(f"[s3-get-object.py] Start getting {object_key}")
    start_time = time.time()

s3 = boto3.client('s3')
with open(outfile, 'wb') as f:
    s3.download_fileobj(BUCKET, object_key, f)

if DEBUG:
    end_time = time.time()
    total_time = end_time - start_time
    total_time = round(total_time, 3)
    print(f"[s3-get-object.py] Finish getting {object_key} in {total_time} seconds")
