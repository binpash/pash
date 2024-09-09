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

def read_graph(filename):
    with open(filename, "rb") as ir_file:
        ir, shell_vars, args = pickle.load(ir_file)
    return ir, shell_vars, args

def invoke_lambda(script_id_to_script, script_id):
    lambda_client = boto3.client("lambda",region_name='us-east-1')
    response = lambda_client.invoke(
        FunctionName="lambda",
        InvocationType="Event",
        LogType="None",
        Payload=json.dumps({"data": json.dumps(script_id_to_script) , "id": script_id}),
    )
    return response

def invoke_lambda_ec2(script_id_to_script, script_id):
    EC2_IP = "129.114.109.71"
    EC2_PORT = 9999
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((EC2_IP, EC2_PORT))
    json_data = json.dumps({"id": script_id, "data": json.dumps(script_id_to_script)})
    encoded_json_data = json_data.encode("utf-8")
    print(len(encoded_json_data))
    # send length first?
    # s.sendall(str(len(encoded_json_data))+"\r\n".encode("utf-8"))
    try:
        s.sendall(json_data.encode("utf-8"))
        s.close()
    except Exception as e:
        print(f"[invoke-lambda.py] EC2 {script_id} invocation error: {e}")
    print(f"[invoke-lambda.py] EC2 {script_id} invocation")

def init(ir_filename: str) -> Tuple[IR, argparse.Namespace, dict]:
    # init pash_args, config, logging
    ir, shell_vars, args = read_graph(ir_filename)
    config.set_config_globals_from_pash_args(args)
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    config.LOGGING_PREFIX = f"Serverless Executor:"
    return ir, shell_vars, args

def main(ir_filename: str):
    ir, shell_vars, args = init(ir_filename)
    # prepare scripts
    main_graph_script_id, main_subgraph_script_id, script_id_to_script = prepare_scripts_for_serverless_exec(ir, shell_vars, args)

    if args.sls_instance == 'lambda':
        response = invoke_lambda(script_id_to_script, main_subgraph_script_id)
        log(response)
        wait_msg_done()
    elif args.sls_instance == 'hybrid':
        invoke_lambda_ec2(script_id_to_script, main_subgraph_script_id)
        wait_msg_done()

if __name__ == '__main__':
    ir_filename = sys.argv[1:][0]
    main(ir_filename)
