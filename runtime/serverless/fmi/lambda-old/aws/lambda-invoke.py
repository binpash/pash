import boto3
import sys
import json

function_name, data, id, num = sys.argv[1:]

lambda_client = boto3.client("lambda")

response = lambda_client.invoke(
    FunctionName=function_name,
    InvocationType="Event",
    LogType="None",
    Payload=json.dumps({"id": id, "data": data, "num": num}),
)
