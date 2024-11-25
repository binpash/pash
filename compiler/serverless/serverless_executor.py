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

def read_graph(filename):
    with open(filename, "rb") as ir_file:
        ir, shell_vars, args = pickle.load(ir_file)
    return ir, shell_vars, args

def init(ir_filename: str) -> Tuple[IR, argparse.Namespace, dict]:
    # init pash_args, config, logging
    ir, shell_vars, args = read_graph(ir_filename)
    config.set_config_globals_from_pash_args(args)
    # config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    # config.LOGGING_PREFIX = f"Serverless Manager:"
    return ir, shell_vars, args

class ThreadSafeCounter:
    def __init__(self, concurrency_limit=1000):
        self.concurrency_limit = concurrency_limit
        self.value = 0
        self.total_time = 0
        self._lock = threading.Lock()
        self.counter_history = []

    def increment(self, num=1):
        # When attempting to increment the counter:
        #   we need to check if the counter is already at the limit
        success = True
        while True:
            if not success:
                time.sleep(1)
            with self._lock:
                if self.value > self.concurrency_limit - num:
                    success = False
                    continue
                self.value += num
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
        # create_dynamo_table()
        self.worker_threads = queue.Queue(maxsize=1000)
        self.lambda_config = BotocoreConfig(
            read_timeout=900,
            connect_timeout=60,
            retries={"max_attempts": 0}
        )
        self.concurrency_limit = 1000
        self.counter = ThreadSafeCounter(self.concurrency_limit)
        self.ec2_enabled = False
        self.ec2_ip = "54.236.13.134"
        self.ec2_port = 9999
    
    def invoke_ec2(self, folder_ids, script_ids):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((self.ec2_ip, self.ec2_port))
        json_data = json.dumps({"ids": script_ids, "folder_ids": folder_ids})
        try:
            s.sendall(json_data.encode("utf-8"))
            s.sendall(b"END")
            _ = s.recv(1024)
            s.close()
            self.counter.decrement()
        except Exception as e:
            log(f"EC2 failed {folder_ids} && {folder_ids}: {e}")

    def invoke_lambda(self, folder_ids, script_ids):
        # log(f"[Serverless Manager] Try to invoke lambda response with batches {script_ids} && {folder_ids}")
        start = time.time()
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
        end = time.time()
        self.counter.decrement(end-start)
        return
    
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
            # prepare scripts
            ir, shell_vars, args = init(ir_filename)
            main_graph_script_id, main_subgraph_script_id, script_id_to_script, ec2_set = prepare_scripts_for_serverless_exec(ir, shell_vars, args, declared_functions_file)
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
            # preempt number of lambdas to be invoked
            self.counter.increment(len(script_id_to_script)-1)
            conn.sendall(bytes_message)
            conn.close()
            for script_id, script in script_id_to_script.items():
                if script_id == main_graph_script_id:
                    continue
                if self.ec2_enabled and script_id in ec2_set:
                    invocation_thread = threading.Thread(target=self.invoke_ec2, args=([s3_folder_id], [script_id]))
                    invocation_thread.start()
                else:
                    invocation_thread = threading.Thread(target=self.invoke_lambda, args=([s3_folder_id], [script_id]))
                    invocation_thread.start()
        
        except Exception as e:
            log(f"[Serverless Manager] Failed to add job to a queue: {e}")
            response_msg = f"ERR {e}"
            bytes_message = response_msg.encode('utf-8')
            conn.sendall(bytes_message)
            conn.close()
