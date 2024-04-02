import boto3
import sys


table_name, primary_key, id, inputfile = sys.argv[1:]

dynamodb = boto3.client("dynamodb")

with open(inputfile, "r") as f:
    data = f.read()

item = {"id": {"S": f"{primary_key}{id}"}, "data": {"S": data}}

try:
    response = dynamodb.put_item(TableName=table_name, Item=item)
except Exception as e:
    print(e)
