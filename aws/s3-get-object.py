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


# import os, sys, time, boto3
# from boto3.s3.transfer import TransferConfig
# from botocore.config import Config

# BUCKET = os.environ["AWS_BUCKET"]
# object_key, outfile = sys.argv[1:]

# # Tune transfer + HTTP pool
# MB = 1024 * 1024
# transfer_cfg = TransferConfig(
#     multipart_threshold=8*MB,        # start multipart at 8MB
#     multipart_chunksize=128*MB,       # part size (try 64–128MB)
#     max_concurrency=32,              # threads per transfer (try 16–32)
#     max_io_queue=2048,               # buffer queue
#     io_chunksize=256*1024,           # per-IO chunk (256KB)
#     num_download_attempts=10
# )

# # Bigger HTTP connection pool for parallel parts
# client = boto3.client(
#     "s3",
#     config=Config(max_pool_connections=128)  # >= 2 * max_concurrency
# )

# with open(outfile, "wb", buffering=0) as f:  # unbuffered to avoid double buffering
#     client.download_fileobj(BUCKET, object_key, f, Config=transfer_cfg)
