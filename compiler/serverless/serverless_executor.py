import argparse
import copy
import json
import os
import pickle
import socket
import sys
import threading
import time
import boto3
from botocore.config import Config as BotocoreConfig
from typing import Any, Tuple
sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
from dspash.socket_utils import SocketManager
from serverless.ir_helper import prepare_scripts_for_serverless_exec
from util import log
from ir import IR
import config
import queue
import multiprocessing

AWS_ACCOUNT_ID=os.environ.get("AWS_ACCOUNT_ID")
QUEUE=os.environ.get("AWS_QUEUE")
EC2_IP=os.environ.get("AWS_EC2_IP")

def wait_msg_done_sqs():
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
        log(f"[Serverless Manager] Table sls-result created!")
    except Exception as e:
        log(f"[Serverless Manager] Table sls-result already exists!")
    return


def wait_msg_done_dynamo(key):
    try:
        while True:
            dynamo = boto3.client('dynamodb')
            table_name = 'sls-result'
            response = dynamo.get_item(
                TableName=table_name,
                Key={
                    'key': {
                        'S': key
                    }
                }
            )
            if 'Item' in response:
                break
            time.sleep(1)    
    except:
        log(f"[Serverless Manager] DynamoDB get item error: {e}")
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
    # config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    # config.LOGGING_PREFIX = f"Serverless Manager:"
    return ir, shell_vars, args

def main(ir_filename: str, declared_functions_file_name: str):

    ir, shell_vars, args = init(ir_filename)
    # prepare scripts
    main_graph_script_id, main_subgraph_script_id, script_id_to_script = prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_file_name)

    # put scripts into s3
    random_id = str(int(time.time()))
    log(f"[Serverless Manager] Uploading scripts to s3 with folder_id: {random_id}")
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
        log(f"[Serverless Manager] Invoke lambda response {response}")
    elif args.sls_instance == 'hybrid':
        invoke_lambda_ec2(script_id_to_script, main_subgraph_script_id, random_id)
    # wait_msg_done()
    wait_msg_done_dynamo(random_id)

def handler(request, conn):
    args = request.split(':', 1)[1].strip()
    ir_filename, declared_functions_file = args.split()

    try:
        ir, shell_vars, args = init(ir_filename)
        # prepare scripts
        main_graph_script_id, main_subgraph_script_id, script_id_to_script = prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_file)

        # put scripts into s3
        random_id = str(int(time.time()))
        s3 = boto3.client('s3')
        bucket = os.getenv("AWS_BUCKET")
        if not bucket:
            raise Exception("AWS_BUCKET environment variable not set")
        for script_id, script in script_id_to_script.items():
            if script_id == main_graph_script_id:
                continue
            s3.put_object(Bucket=bucket, Key=f'sls-scripts/{random_id}/{script_id}.sh', Body=script)

        if args.sls_instance == 'lambda':
            response = invoke_lambda(script_id_to_script, main_subgraph_script_id, random_id)
            log(f"[Serverless Manager] Invoke lambda response {response}")
        elif args.sls_instance == 'hybrid':
            invoke_lambda_ec2(script_id_to_script, main_subgraph_script_id, random_id)

        response_msg = f"OK {random_id}"
        bytes_message = response_msg.encode('utf-8')
        conn.sendall(bytes_message)
        conn.close()
    
    except Exception as e:
        log(f"[Serverless Manager] Executing lambda error {e}")
        response_msg = f"ERR {e}"
        bytes_message = response_msg.encode('utf-8')
        conn.sendall(bytes_message)
        conn.close()
class ThreadSafeCounter:
    def __init__(self, concurrency_limit=1000):
        self.concurrency_limit = concurrency_limit
        self.value = 0
        self.total_time = 0
        self._lock = threading.Lock()
        self.counter_history = []

    def increment(self):
        # When attempting to increment the counter:
        #   we need to check if the counter is already at the limit
        success = True
        while True:
            if not success:
                time.sleep(1)
            with self._lock:
                if self.value == self.concurrency_limit:
                    success = False
                    continue
                self.value += 1
                # log(f"[Serverless Manager] Counter value: {self.value}")
                success = True
                self.counter_history.append((time.time(), self.value))
                break
        return self.value
    
    def decrement(self, time_taken=0):
        with self._lock:
            self.value -= 1
            self.total_time += time_taken
            self.counter_history.append((time.time(), self.value))
        return self.value

    def get_value(self):
        with self._lock:
            return self.value, self.total_time
    
    def print_history(self):
        with self._lock:
            for time, value in self.counter_history:
                log(f"[Serverless Manager] Counter value: {value} at time: {time}")
        return

class ServerlessManager:
    def __init__(self):
        create_dynamo_table()
        self.worker_threads = queue.Queue(maxsize=1000)
        self.lambda_config = BotocoreConfig(
            read_timeout=900,
            connect_timeout=60,
            retries={"max_attempts": 0}
        )
        self.concurrency_limit = 1000
        self.counter = ThreadSafeCounter(self.concurrency_limit)

    def invoke_lambda(self, folder_ids, script_ids):
        # log(f"[Serverless Manager] Try to invoke lambda response with batches {script_ids} && {folder_ids}")
        try:
            lambda_client = boto3.client("lambda",region_name='us-east-1', config=self.lambda_config)
            response = lambda_client.invoke(
                FunctionName="lambda",
                InvocationType="RequestResponse",
                LogType="None",
                Payload=json.dumps({"folder_ids": folder_ids , "ids": script_ids}),
            )
            # parse response
            # log(f"[Serverless Manager] Lambda execution response with batches {script_ids} && {folder_ids}: {response}")
            if response['StatusCode'] != 200:
                log(f"[Serverless Manager] Lambda execution failed with batches {script_ids} && {folder_ids}: {response}")
        except Exception as e:
            log(f"[Serverless Manager] Lambda execution raises exception with batches {script_ids} && {folder_ids}: {e}")
        return
    
    def force(self):
        folder_ids, script_ids = [], []
        last_time_forced = time.time()
        while True and not self.exit:
            if self.job_queue.empty():
                time.sleep(1)
                continue
            s3_folder_id, script_id = self.job_queue.get()
            folder_ids.append(s3_folder_id)
            script_ids.append(script_id)
            # force the lambda with either a batch or after one sec
            if len(folder_ids) == self.job_batch_size or time.time() - last_time_forced > 2:
                random_id = str(int(time.time()))
                worker_thread = threading.Thread(target=self.invoke_lambda, args=(folder_ids, script_ids))
                worker_thread.start()
                self.worker_threads.put(worker_thread)
                folder_ids, script_ids = [], []
                last_time_forced = time.time()

    def run(self):
        self.sls_socket = SocketManager(os.getenv('SERVERLESS_SOCKET'))

        while True:
            request, conn = self.sls_socket.get_next_cmd()

            if request.startswith("Done"):
                counter, total_time = self.counter.get_value()
                while counter > 0:
                    time.sleep(1)
                    counter, total_time = self.counter.get_value()
                self.counter.print_history()
                log(f"[Serverless Manager] All jobs are done with total time: {total_time}")
                response_msg = f"All jobs are done, shutting down the serverless manager"
                self.sls_socket.respond(response_msg, conn)
                self.sls_socket.close()
                break

            # multi-threaded
            client_thread = threading.Thread(target=self.handler, args=(request, conn))
            client_thread.start()

    def handler(self, request, conn):
        # print(f"[Serverless Manager] Received request: {request}")
        args = request.split(':', 1)[1].strip()
        ir_filename, declared_functions_file = args.split()

        try:
            now = time.time()
            ir, shell_vars, args = init(ir_filename)
            # prepare scripts
            main_graph_script_id, main_subgraph_script_id, script_id_to_script = prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_file)

            # put scripts into s3
            s3_folder_id = str(int(time.time()))
            s3 = boto3.client('s3')
            bucket = os.getenv("AWS_BUCKET")
            if not bucket:
                raise Exception("AWS_BUCKET environment variable not set")
            for script_id, script in script_id_to_script.items():
                if script_id == main_graph_script_id:
                    continue
                s3.put_object(Bucket=bucket, Key=f'sls-scripts/{s3_folder_id}/{script_id}.sh', Body=script)

            response_msg = f"OK {s3_folder_id} {main_subgraph_script_id}"
            bytes_message = response_msg.encode('utf-8')
            conn.sendall(bytes_message)
            conn.close()
            start = time.time()
            self.counter.increment()
            self.invoke_lambda([s3_folder_id], [main_subgraph_script_id])
            end = time.time()
            if len(script_id_to_script) > 2:
                wait_msg_done_dynamo(s3_folder_id)
            self.counter.decrement(time_taken=end-start)
        
        except Exception as e:
            log(f"[Serverless Manager] Failed to add job to a queue: {e}")
            response_msg = f"ERR {e}"
            bytes_message = response_msg.encode('utf-8')
            conn.sendall(bytes_message)
            conn.close()

if __name__ == '__main__':
    ir_filename = sys.argv[1:][0]
    declared_functions =  sys.argv[1:][1]
    main(ir_filename, declared_functions)
