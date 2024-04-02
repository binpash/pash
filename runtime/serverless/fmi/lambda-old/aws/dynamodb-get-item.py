import boto3
import sys

table_name, primary_key, id, outfile = sys.argv[1:]

dynamodb = boto3.client("dynamodb")

try:
    response = dynamodb.get_item(
        TableName=table_name, Key={"id": {"S": f"{primary_key}{id}"}}
    )

    item = response["Item"]

    with open(outfile, "w") as f:
        print(item["data"]["S"], file=f, end="")

except Exception as e:
    print(e)
