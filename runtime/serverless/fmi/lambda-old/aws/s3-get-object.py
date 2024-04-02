import boto3
import sys

bucket_name, object_key, outfile = sys.argv[1:]

session = boto3.Session()

s3 = session.client("s3")

response = s3.get_object(Bucket=bucket_name, Key=object_key)

with open(outfile, "w") as f:
    print(response["Body"].read().decode("utf-8"), file=f, end="")
