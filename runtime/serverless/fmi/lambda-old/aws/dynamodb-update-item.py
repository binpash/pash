import boto3
import sys


table_name, primary_key, id = sys.argv[1:]

dynamodb = boto3.client("dynamodb")

update_expression = f"SET #attr = :val"

expression_attribute_names = {"#attr": f"w{id}"}

expression_attribute_values = {":val": {"S": "1"}}

try:
    response = dynamodb.update_item(
        TableName=table_name,
        Key={"id": {"S": primary_key}},
        UpdateExpression=update_expression,
        ExpressionAttributeNames=expression_attribute_names,
        ExpressionAttributeValues=expression_attribute_values,
        ReturnValues="ALL_NEW",
    )
except Exception as e:
    print(e)

should_aggregate = True

for k, v in response["Attributes"].items():
    if k == "id":
        continue
    if v["S"] == "0":
        should_aggregate = False
        break

if should_aggregate:
    print(1)
else:
    print(0)
