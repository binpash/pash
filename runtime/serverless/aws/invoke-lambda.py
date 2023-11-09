import boto3
import sys
import json
print(sys.argv[1:])
id_, data_path = sys.argv[1:]
with open(data_path, "r") as f:
    data = f.read()

lambda_client = boto3.client("lambda")

response = lambda_client.invoke(
    FunctionName="lambda",
    InvocationType="Event",
    LogType="None",
    Payload=json.dumps({"id": id_, "data": data}),
)
