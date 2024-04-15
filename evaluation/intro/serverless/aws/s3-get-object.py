import boto3
import sys
import time
import os

BUCKET=os.environ.get("AWS_BUCKET")
object_key, outfile = sys.argv[1:]

print(f"[s3-get-object.py] Start getting {object_key}")
start_time = time.time()

session = boto3.Session()
s3 = session.client("s3")

try:
    response = s3.get_object(Bucket=BUCKET, Key=object_key)
except Exception as e:
    print(e)

batch = 10000

with open(outfile, "wb") as f:
    while True:
        x = response["Body"].read(batch)

        if not x:
            break

        f.write(x)
        f.flush()

end_time = time.time()
total_time = end_time - start_time
total_time = round(total_time, 3)

print(f"[s3-get-object.py] Finish getting {object_key} in {total_time} seconds")
