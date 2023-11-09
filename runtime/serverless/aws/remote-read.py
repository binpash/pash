import boto3
import sys
import time

object_key, outfile = sys.argv[1:]
print("Remote read", object_key)
session = boto3.Session()

s3 = session.client("s3")

while True:
    try:
        response = s3.get_object(Bucket="yizhengx", Key=object_key)
        break
    except:
        time.sleep(1)

with open(outfile, "wb") as f:
    f.write(response["Body"].read())