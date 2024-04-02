import boto3
import sys

bucket_name, object_key, infile = sys.argv[1:]

session = boto3.Session()

s3 = session.client("s3")

# Prepare the data to be stored in the S3 object
with open(infile, "rb") as file:
    object_data = file.read()

# Upload the object to S3
s3.put_object(Bucket=bucket_name, Key=object_key, Body=object_data)
