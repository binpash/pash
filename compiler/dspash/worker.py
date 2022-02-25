import socket, pickle
from threading import Thread
from socket_utils import encode_request, decode_request
import subprocess
import json
import sys
import os
import asyncio
import shlex
import time
import uuid
import argparse

PASH_TOP = os.environ['PASH_TOP']
sys.path.append(os.path.join(PASH_TOP, "compiler"))

import config
from util import log
from annotations import load_annotation_files
import pash_runtime
from dspash.socket_utils import send_msg, recv_msg
from dspash.ir_helper import graph_to_shell

# from ... import config
HOST = '0.0.0.0'
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)


EXEC_SCRIPT_PATH = f'{PASH_TOP}/compiler/dspash/remote_exec_script.sh'

def err_print(*args):
    print(*args, file=sys.stderr)

def send_success(conn, body, msg = ""):
    request = {
        'status': 'OK',
        'body': body,
        'msg': msg
    }
    send_msg(conn, encode_request(request))

def parse_exec_request(request):
    return request['cmd']

def parse_exec_graph(request):
    return request['graph'], request['shell_variables']

def exec_graph(graph, shell_vars):
    config.config['shell_variables'] = shell_vars
    script_path = graph_to_shell(graph)

    e = os.environ.copy()
    e['PASH_TOP'] = PASH_TOP
    rc = subprocess.Popen(['bash', script_path], env=e)
    return rc

class Worker:
    def __init__(self, port = None):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if port == None:
            # pick a random port
            self.s.bind((HOST, 0))
        else:
            self.s.bind((HOST, port))
        print(f"Worker running on port {self.s.getsockname()[1]}")

    def run(self):
        connections = []
        with self.s:
            self.s.listen()
            while(True):
                conn, addr = self.s.accept()
                print(f"got new connection")     
                t = Thread(target=manage_connection, args=[conn, addr])
                t.start()
                connections.append(t)
        for t in connections:
            t.join()

def manage_connection(conn, addr):
    rcs = []
    with conn:
        print('Connected by', addr)
        while True:
            data = recv_msg(conn)
            if not data:
                break

            print("got new request")
            request = decode_request(data)      
            # if request['type'] == 'exec':
            #     cmd = parse_exec_request(request)
            #     from_port = get_available_port()
            #     to_port = get_available_port()
            #     print(cmd, from_port, to_port)
            #     args = ["bash", EXEC_SCRIPT_PATH, addr[0], from_port, to_port, cmd]
            #     rc = subprocess.Popen(args)
            #     body = {
            #         'from_port': from_port,
            #         'to_port': to_port
            #     }
            if request['type'] == 'Exec-Graph':
                graph, shell_vars = parse_exec_graph(request)
                exec_graph(graph, shell_vars)
                body = {}
            else:
                print(f"Unsupported request {request}")
            send_success(conn, body)
    print("connection ended")
    for rc in rcs:
        rc.wait()

def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--port",
                        type=int,
                        help="port to use",
                        default=65432)
    config.add_common_arguments(parser)
    args = parser.parse_args()
    config.pash_args = args
    ## Initialize the log file
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    return args

def init():
    args = parse_args()
    config.LOGGING_PREFIX = f"Worker {config.pash_args.port}: "
    config.annotations = load_annotation_files(
        config.config['distr_planner']['annotations_dir'])
    pash_runtime.runtime_config = config.config['distr_planner']

def main():
    init()
    worker = Worker(config.pash_args.port)
    worker.run()

if __name__ == "__main__":
    main()
    
