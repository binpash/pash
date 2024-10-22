import boto3
import sys
import json
id_, data_ = sys.argv[1:]
print("[invoke-lambda.py] Invoke lambda", id_)

lambda_client = boto3.client("lambda")

response = lambda_client.invoke(
    FunctionName="lambda",
    InvocationType="Event",
    LogType="None",
    Payload=json.dumps({"id": id_, "data": data_}),
)
