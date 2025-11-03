import boto3
import sys
import json
import socket

id_, folder_id, type_ = sys.argv[1:]
EC2_IP = "34.228.68.128"
EC2_PORT = 9999

if type_ == "lambda":
    print(f"[invoke-lambda.py] Try invoking lambda {id_}")
    lambda_client = boto3.client("lambda")
    response = lambda_client.invoke(
        FunctionName="lambda",
        InvocationType="Event",
        LogType="None",
        Payload=json.dumps({"ids": [id_], "folder_ids": [folder_id]}),
    )
    print(f"[invoke-lambda.py] Invoked lambda {id_}: {response}")
elif type_ == "ec2":
    print(f"[invoke-lambda.py] Invoke ec2 {id_}")
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((EC2_IP, EC2_PORT))
    try:
        json_str = json.dumps({"id": id_, "folder_id": folder_id})
        s.sendall(json_str.encode("utf-8"))
        # response = s.recv(1024)
        s.close()
        print(f"[invoke-lambda.py] EC2 {id_} invocation sent")
    except Exception as e:
        print(f"[invoke-lambda.py] EC2 {id_} invocation error: {e}")