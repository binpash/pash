import signal
import socket
from multiprocessing import Process
import time
from socket_utils import encode_request, decode_request
import subprocess
import sys
import os
import argparse
import requests

DISH_TOP = os.environ['DISH_TOP']
PASH_TOP = os.environ['PASH_TOP']
sys.path.append(os.path.join(PASH_TOP, "compiler"))

from dspash.utils import create_filename, write_file
from dspash.ir_helper import save_configs, to_shell_file, to_shell, add_debug_flags, add_kill_flags
from dspash.socket_utils import send_msg, recv_msg
import pash_runtime
from annotations import load_annotation_files
from util import log
import config

# from ... import config
HOST = socket.gethostbyname(socket.gethostname())
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)


def err_print(*args):
    print(*args, file=sys.stderr, flush=True)


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


def exec_graph(graph, shell_vars, functions, kill_target, debug=False):
    config.config['shell_variables'] = shell_vars
    if debug:
        err_print('debug is on')
        add_debug_flags(graph)
        stderr = subprocess.PIPE
    else:
        stderr = None

    if kill_target:
        log(f'going to kill: {kill_target}')
        add_kill_flags(graph, kill_target)

    script_path = to_shell_file(graph, config.pash_args)

    e = os.environ.copy()
    e['PASH_TOP'] = PASH_TOP
    e['DISH_TOP'] = DISH_TOP

    # store functions
    functions_file = create_filename(
        dir=config.PASH_TMP_PREFIX, prefix='pashFuncs')
    write_file(functions_file, functions)
    cmd = f"source {functions_file}; source {script_path}"
    err_print(f"executing {cmd}")

    rc = subprocess.Popen(
        cmd, env=e, executable="/bin/bash", shell=True, stderr=stderr, preexec_fn=os.setsid)
    return rc


class Worker:
    def __init__(self, port=None):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        for _ in range(5):
            try:
                if port == None:
                    # pick a random port
                    self.s.bind((HOST, 0))
                else:
                    self.s.bind((HOST, port))
                break
            except socket.error as e:
                err_print(f"Bind failed with error: {e}. Retrying...")
                time.sleep(3)
        else:
            raise Exception("Could not bind after multiple attempts")
        err_print(f"Worker running on {HOST}:{self.s.getsockname()[1]}")

    def run(self):
        connections = []
        with self.s:
            self.s.listen()
            while (True):
                conn, addr = self.s.accept()
                err_print(f"got new connection")
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
        err_print("process timedout")
        rc.kill()
        _, err = rc.communicate()

    response = {
        'name': name,
        'returncode': rc.returncode,
        'stderr': err.decode("UTF-8"),
        'shellscript': shell_script,
    }

    requests.post(url=url, json=response)


def manage_connection(conn, addr):
    rcs = []
    with conn:
        err_print('Connected by', addr)
        dfs_configs_paths = {}
        while True:
            data = recv_msg(conn)
            if not data:
                break

            request = decode_request(data)
            err_print("got new request", request['type'])

            if request['type'] == 'Exec-Graph':
                graph, shell_vars, functions = parse_exec_graph(request)
                debug = True if request['debug'] else False
                kill_target = request['kill_target']
                save_configs(graph, dfs_configs_paths)
                rc = exec_graph(graph, shell_vars, functions, kill_target, debug)
                rcs.append((rc, request))
                body = {}
            elif request['type'] == 'Kill-All-Subgraphs':
                body = handle_kill_all_subgraphs_request(rcs)
            else:
                err_print(f"Unsupported request {request}")
            send_success(conn, body)
    err_print("connection ended")

    err_print("len rcs:", len(rcs))
    for rc, request in rcs:
        err_print("rc")
        if request['debug']:
            send_log(rc, request)
        else:
            rc.wait()


def handle_kill_all_subgraphs_request(rcs):
    # Record the start time
    start_time = time.time()
    err_print("Killing every subgraph")

    # Try to terminate the subprocess gracefully
    for rc, _ in rcs:
        os.killpg(os.getpgid(rc.pid), signal.SIGTERM)
        # rc.terminate()

    for rc, _ in rcs:
        try:
            # Wait for the subprocess to terminate, with a timeout
            rc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            # If the subprocess is still running after the timeout, kill it
            err_print(f"Terminating subprocess {rc.pid} failed, killing it forcefully")
            os.killpg(os.getpgid(rc.pid), signal.SIGKILL)
            rc.wait()

    # Record the end time
    end_time = time.time()

    # Calculate the time elapsed
    elapsed_time = end_time - start_time

    err_print(f"Killed every subgraph in {elapsed_time} seconds")
    return {}


def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--port",
                        type=int,
                        help="port to use",
                        default=65432)
    config.add_common_arguments(parser)
    args = parser.parse_args()
    config.pash_args = args
    # Initialize the log file
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
    # needed because graphs could have multiples sinks
    config.pash_args.termination = ""


def main():
    init()
    worker = Worker(config.pash_args.port)
    worker.run()


if __name__ == "__main__":
    main()
