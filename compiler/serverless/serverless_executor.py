import argparse
import copy
import json
import os
import pickle
import socket
import sys
import time
import boto3
from typing import Any, Tuple
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
from serverless.ir_helper import prepare_scripts_for_serverless_exec
from util import log
from ir import IR
import config

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
QUEUE=os.environ.get("AWS_QUEUE")
EC2_IP=os.environ.get("AWS_EC2_IP")

def exec():
    pass

def wait_msg_done():
    while True:
        sqs = boto3.client('sqs')
        queue_url = f'https://sqs.us-east-1.amazonaws.com/{AWS_ACCOUNT_ID}/{QUEUE}'
        # Receive message from SQS queue
        response = sqs.receive_message(
            QueueUrl=queue_url,
            AttributeNames=[
                'SentTimestamp'
            ],
            MaxNumberOfMessages=1,
            MessageAttributeNames=[
                'All'
            ],
            VisibilityTimeout=30,
            WaitTimeSeconds=20
        )
        try:
            message = response['Messages'][0]
            receipt_handle = message['ReceiptHandle']

            # Delete received message from queue
            sqs.delete_message(
                QueueUrl=queue_url,
                ReceiptHandle=receipt_handle
            )
            log('Received and deleted message: %s' % message)
            break
        except:
            time.sleep(1)

def create_dynamo_table():
    dynamo = boto3.client('dynamodb')
    table_name = 'sls-result'
    try:
        response = dynamo.create_table(
            TableName=table_name,
            KeySchema=[
                {
                    'AttributeName': 'key',
                    'KeyType': 'HASH'
                }
            ],
            AttributeDefinitions=[
                {
                    'AttributeName': 'key',
                    'AttributeType': 'S'
                }
            ],
            BillingMode='PAY_PER_REQUEST'
        )
        log(f"Table sls-result created!")
    except Exception as e:
        log(f"Table sls-result already exists!")
    return


def wait_msg_done_dynamo(key):
    while True:
        dynamo = boto3.client('dynamodb')
        table_name = 'sls-result'
        try:
            response = dynamo.get_item(
                TableName=table_name,
                Key={
                    'key': {
                        'S': key
                    }
                }
            )
            item = response['Item']
            log('Received message: %s' % item)
            break
        except:
            # log("Waiting for message")
            time.sleep(1)
    return

def write_to_dynamo(key, value):
    dynamo = boto3.client('dynamodb')
    table_name = 'sls-result'
    try:
        response = dynamo.put_item(
            TableName=table_name,
            Item={
                'key': {
                    'S': key
                },
                'output_id': {
                    'S': value
                }
            }
        )
        log(f"Successfully wrote item {key}:{value} to dynamoDB.")
    except Exception as e:
        log(f"Failed to write item {key}:{value} to dynamoDB.")
    return

def read_graph(filename):
    with open(filename, "rb") as ir_file:
        ir, shell_vars, args = pickle.load(ir_file)
    return ir, shell_vars, args

def invoke_lambda(script_id_to_script, script_id, folder_id):
    lambda_client = boto3.client("lambda",region_name='us-east-1')
    response = lambda_client.invoke(
        FunctionName="lambda",
        InvocationType="Event",
        LogType="None",
        Payload=json.dumps({"folder_id": folder_id , "id": script_id}),
    )
    return response

def invoke_lambda_ec2(script_id_to_script, script_id, folder_id):
    EC2_PORT = 9999
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((EC2_IP, EC2_PORT))
    json_data = json.dumps({"id": script_id, "folder_id": folder_id})
    # send length first?
    # s.sendall(str(len(encoded_json_data))+"\r\n".encode("utf-8"))
    try:
        s.sendall(json_data.encode("utf-8"))
        s.close()
    except Exception as e:
        log(f"EC2 {script_id} invocation error: {e}")
    log(f"EC2 {script_id} invocation")

def init(ir_filename: str) -> Tuple[IR, argparse.Namespace, dict]:
    # init pash_args, config, logging
    ir, shell_vars, args = read_graph(ir_filename)
    config.set_config_globals_from_pash_args(args)
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    config.LOGGING_PREFIX = f"Serverless Executor:"
    return ir, shell_vars, args

def main(ir_filename: str, declared_functions_file_name: str):

    ir, shell_vars, args = init(ir_filename)
    # prepare scripts
    main_graph_script_id, main_subgraph_script_id, script_id_to_script = prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_file_name)

    # put scripts into s3
    random_id = str(int(time.time()))
    log(f"Uploading scripts to s3 with folder_id: {random_id}")
    s3 = boto3.client('s3')
    bucket = os.getenv("AWS_BUCKET")
    if not bucket:
        raise Exception("AWS_BUCKET environment variable not set")
    for script_id, script in script_id_to_script.items():
        s3.put_object(Bucket=bucket, Key=f'sls-scripts/{random_id}/{script_id}.sh', Body=script)

    # create dynamo table if not exists
    create_dynamo_table()

    if args.sls_instance == 'lambda':
        response = invoke_lambda(script_id_to_script, main_subgraph_script_id, random_id)
        log(response)
    elif args.sls_instance == 'hybrid':
        invoke_lambda_ec2(script_id_to_script, main_subgraph_script_id, random_id)
    # wait_msg_done()
    wait_msg_done_dynamo(random_id)

if __name__ == '__main__':
    ir_filename = sys.argv[1:][0]
    declared_functions =  sys.argv[1:][1]
    main(ir_filename, declared_functions)
