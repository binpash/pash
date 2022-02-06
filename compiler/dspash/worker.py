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
# from ... import config
HOST = '0.0.0.0'
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)

#TODO: figure how to import from top directory or move files up
if 'PASH_TOP' in os.environ:
    PASH_TOP = os.environ['PASH_TOP']
else:
    GIT_TOP_CMD = [ 'git', 'rev-parse', '--show-toplevel', '--show-superproject-working-tree']
    PASH_TOP = subprocess.run(GIT_TOP_CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.rstrip()

EXEC_SCRIPT_PATH = f'{PASH_TOP}/compiler/dspash/remote_exec_script.sh'

def err_print(*args):
    print(*args, file=sys.stderr)

def send_success(conn, body, msg = ""):
    request = {
        'status': 'OK',
        'body': body,
        'msg': msg
    }
    conn.send(encode_request(request))

def parse_exec_request(request):
    return request['cmd']


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


def get_available_port():
    # There is a possible race condition using the returned port as it could be opened by a different process
    cmd = "bash -c \"comm -23 <(seq 1100 65535 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | shuf | head -n 1\""
    args = shlex.split(cmd)
    p = subprocess.run(args, stdout=subprocess.PIPE, universal_newlines=True)
    return p.stdout.rstrip()


def manage_connection(conn, addr):
    rcs = []
    with conn:
        print('Connected by', addr)
        while True:
            request = conn.recv(4096)
            if not request:
                break

            print("got new request")
            request = decode_request(request)
            if request['type'] == 'exec':
                cmd = parse_exec_request(request)
                from_port = get_available_port()
                to_port = get_available_port()
                print(cmd, from_port, to_port)
                args = ["bash", EXEC_SCRIPT_PATH, addr[0], from_port, to_port, cmd]
                e = os.environ.copy()
                e['PASH_TOP'] = PASH_TOP
                rc = subprocess.Popen(args, env = e)
                body = {
                    'from_port': from_port,
                    'to_port': to_port
                }
                time.sleep(0.1) # kinda hacky but sometimes port is not open by time host gets response
                send_success(conn, body)
    print("connection ended")
    for rc in rcs:
        rc.wait()

def main():
    connections = []
    port = PORT if len(sys.argv) < 2 else int(sys.argv[1])
    worker = Worker(port)
    worker.run()

if __name__ == "__main__":
    main()
    
