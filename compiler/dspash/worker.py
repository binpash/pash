import socket
from multiprocessing import Process
from socket_utils import encode_request, decode_request
import subprocess
import sys
import os
import argparse
import requests
import time
import threading
import signal

DISH_TOP = os.environ['DISH_TOP']
PASH_TOP = os.environ['PASH_TOP']
sys.path.append(os.path.join(PASH_TOP, "compiler"))

from dspash.utils import create_filename, write_file
from dspash.ir_helper import save_configs, to_shell_file, to_shell, add_debug_flags
from dspash.socket_utils import send_msg, recv_msg
import pash_compiler
from util import log
import config

# from ... import config
HOST = socket.gethostbyname(socket.gethostname())
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)

def err_print(*args):
    print(*args, file=sys.stderr)


def send_success(conn, body, msg=""):
    request = {
        'status': 'OK',
        'body': body,
        'msg': msg
    }
    send_msg(conn, encode_request(request))


def parse_exec_request(request):
    return request['cmd']


def parse_exec_graph(request):
    return request['graph'], request['shell_variables'], request['functions']


def exec_graph(graph, shell_vars, functions, debug=False):
    config.config['shell_variables'] = shell_vars
    if debug:
        log('debug is on')
        add_debug_flags(graph)
        stderr = subprocess.PIPE
    else:
        stderr = None

    script_path = to_shell_file(graph, config.pash_args)
    e = os.environ.copy()
    e['PASH_TOP'] = PASH_TOP
    e['DISH_TOP'] = DISH_TOP
    # store functions
    functions_file = create_filename(
        dir=config.PASH_TMP_PREFIX, prefix='pashFuncs')
    write_file(functions_file, functions)
    cmd = f"source {functions_file}; source {script_path}"

    log(f"Executing: {cmd}")
    rc = subprocess.Popen(
        cmd, env=e, executable="/bin/bash", shell=True, stderr=stderr, preexec_fn=os.setsid)
    return rc


class Worker:
    def __init__(self, port=None):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if port == None:
            # pick a random port
            self.s.bind((HOST, 0))
        else:
            self.s.bind((HOST, port))
        log(f"Worker running on {HOST}:{self.s.getsockname()[1]}")

    def run(self):
        connections = []
        with self.s:
            self.s.listen()
            while (True):
                conn, addr = self.s.accept()
                log(f"got new connection")
                t = Process(target=manage_connection, args=[conn, addr])
                t.start()

                connections.append(t)
        for t in connections:
            t.join()


def send_log(rc: subprocess.Popen, request):
    name = request['debug']['name']
    url = request['debug']['url']
    shell_script = to_shell(request['graph'], config.pash_args)

    try:
        # timeout is set to 10s for debuggin
        _, err = rc.communicate(timeout=10)
    except:
        log("process timedout")
        rc.kill()
        _, err = rc.communicate()

    response = {
        'name': name,
        'returncode': rc.returncode,
        'stderr': err.decode("UTF-8"),
        'shellscript': shell_script,
    }

    requests.post(url=url, json=response)

def send_discovery_server_log(rc: subprocess.Popen, request):
    name = f"{request['debug']['name']}:Discovery Server"
    url = request['debug']['url']
    print(type(rc), rc)
    log_output = rc.stderr.read()

    response = {
        'name': name,
        'returncode': 0,
        'stderr': log_output.decode("UTF-8") if log_output else "N/A",
        'shellscript': 'N/A'
    }

    requests.post(url=url, json=response)


def manage_connection(conn, addr):
    rcs = []
    with conn:
        log('Connected by', addr)
        dfs_configs_paths = {}
        while True:
            try:
                data = recv_msg(conn)
                if not data:
                    break
                log("got new request")
                request = decode_request(data) 
                log(request)           
                if request['type'] == 'Exec-Graph':
                    graph, shell_vars, functions = parse_exec_graph(request)
                    debug = True if request['debug'] else False
                    save_configs(graph, dfs_configs_paths)
                    time.sleep(int(request['worker_timeout']))
                    rc = exec_graph(graph, shell_vars, functions, debug)
                    rcs.append((rc, request))
                    body = {}
                elif request['type'] == 'Done':
                    log("Received 'Done' signal. Closing connection from the worker.")
                    break
                elif request['type'] == 'abortAll':
                    # This is buggy so not used
                    num_aborted = 0
                    while rcs:
                        rc, request = rcs.pop()
                        os.killpg(os.getpgid(rc.pid), signal.SIGTERM)

                        # os.system('pkill -TERM -P {pid}'.format(pid=rc.pid))
                        num_aborted += 1
                    body = {"num_aborted": num_aborted}
                else:
                    log(f"Unsupported request {request}")
                send_success(conn, body)
            except Exception as e:
                log(e)
                break

    # Ensure subprocesses have finished, and releasing corresponding resources

    for rc, request in rcs:
        if request['debug']:
            send_log(rc, request)
        else:
            rc.wait()
    log("connection ended")


def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--port",
                        type=int,
                        help="port to use",
                        default=65432)
    config.add_common_arguments(parser)
    args = parser.parse_args()
    config.set_config_globals_from_pash_args(args)
    ## Initialize the log file
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)
    return args


def init():
    args = parse_args()
    config.LOGGING_PREFIX = f"Worker {config.pash_args.port}: "
    ## KK: 2023-02-21 Commenting this out, we need to figure out if the new annotations work with the distribution package
    # config.annotations = load_annotation_files(
    #     config.config['distr_planner']['annotations_dir'])
    pash_compiler.runtime_config = config.config['distr_planner']
    pash_compiler.termination = ""


def main():
    init()
    # start discovery server and track its logs
    # dish_top = os.getenv('DISH_TOP')
    # command = [f'{dish_top}/runtime/dspash/file_reader/discovery_server', '&']
    # discovery_server = subprocess.Popen(command, stderr=subprocess.PIPE, shell=True)

    worker = Worker(config.pash_args.port)
    worker.run()


if __name__ == "__main__":
    main()
