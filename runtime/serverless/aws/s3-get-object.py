import boto3
import sys
import time

object_key, outfile = sys.argv[1:]
print("[s3-get-object.py] Remote read", object_key)
session = boto3.Session()

s3 = session.client("s3")

try:
    response = s3.get_object(Bucket="bucket", Key=object_key)
except Exception as e:
    print(e)

batch = 1000

with open(outfile, "wb") as f:
    while True:
        x = response["Body"].read(batch)

        if not x:
            break

        f.write(x)
        f.flush()
