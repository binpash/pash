import argparse
import copy
import json
import os
import pickle
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

    # start serverless execution by invoking the first lambda
    response = invoke_lambda(script_id_to_script, main_subgraph_script_id)
    log(response)
    wait_msg_done()

if __name__ == '__main__':
    ir_filename= args = sys.argv[1:][0]
    main(ir_filename)
