import shutil
import signal
import socket
import time
from typing import List
from socket_utils import encode_request, decode_request
import subprocess
import sys
import os
import argparse
import requests
import time
import stat
from threading import Event, Thread
from collections import deque

DISH_TOP = os.environ['DISH_TOP'] 
PASH_TOP = os.environ['PASH_TOP']
sys.path.append(os.path.join(PASH_TOP, "compiler"))

from dspash.utils import create_filename, write_file
from dspash.ir_helper import save_configs, to_shell_file, to_shell, add_debug_flags, add_kill_flags
from dspash.socket_utils import send_msg, recv_msg
import pash_runtime
from annotations import load_annotation_files
import config


HOST = socket.gethostbyname(socket.gethostname())
PORT = 65432
DEBUG = True


class Worker:
    def __init__(self):
        self.worker_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.worker_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.last_exec_time_dict = {"": 5000}
        for _ in range(5):
            try:
                self.worker_socket.bind((HOST, PORT))
                err_print(f"Worker running on {HOST}:{self.worker_socket.getsockname()[1]}")
                break
            except socket.error as e:
                err_print(f"Bind failed with error: {e}. Retrying...")
                time.sleep(3)
        else:
            raise Exception("Could not bind after multiple attempts")

    def run(self):
        with self.worker_socket:
            self.worker_socket.listen()
            cnt = 1
            while (True):
                conn, addr = self.worker_socket.accept()
                name = "RequestHandler-" + str(cnt)
                cnt += 1
                err_print("Got new connection running", name)
                RequestHandler(name, conn, addr, self).start()


class RequestHandler(Thread):
    def __init__(self, name, conn, addr, worker: Worker):
        Thread.__init__(self)
        self.name = name
        self.conn = conn
        self.addr = addr
        self.rc_graph_merger_list = []
        self.dfs_configs_paths = {}
        self.body = {}
        self.worker = worker
        self.env_var = os.environ.copy()
        self.env_var['PASH_TOP'] = PASH_TOP
        self.env_var['DISH_TOP'] = DISH_TOP
        self.ft = ""
        self.debug = 0
        self.kill_target = None
        self.event_loop = None
        self.first_request = True

    def run(self):
        with self.conn:
            err_print('Connected by', self.addr, 'name', self.name)
            start_time = time.time()
            doCleanup = False

            while True:
                data = recv_msg(self.conn)
                if not data:
                    break
                self.request = decode_request(data)
                err_print("Got a new request", self.request['type'])

                if self.request['type'] == 'Setup':
                    self.handle_setup_request()
                    doCleanup = True
                elif self.request['type'] == 'Kill-Node':
                    self.handle_kill_node()
                    doCleanup = True
                elif self.request['type'] == 'Exec-Graph':
                    self.handle_exec_graph_request()
                    doCleanup = True
                elif self.request['type'] == 'Batch-Exec-Graph':
                    self.handle_batch_exec_graph()
                    doCleanup = True
                elif self.request['type'] == 'Kill-Subgraphs':
                    self.handle_kill_subgraphs_request()
                    doCleanup = True
                elif self.request['type'] == 'resurrect':
                    self.handle_resurrect_request()
                    # No need to clean up if we just want to bring up the current node
                    doCleanup = False
                else:
                    err_print(f"Unsupported request {self.request}")
                send_success(self.conn, self.body)
            if doCleanup:
                self.cleanup(start_time)

    def cleanup(self, start_time):
        err_print("Connection ended. Subprocess count:", len(self.rc_graph_merger_list))

        for rc, _, _ in self.rc_graph_merger_list:
            rc.wait()

        # record time before sending logs
        end_time = time.time()
        if self.ft == "optimized":
            self.event_loop.quit.set()
            self.event_loop.join()
            err_print(f"{self.event_loop.name} joined")

        # Skip this for now, we no longer send logs to the flask app
        # if self.debug:
        #     for rc, script_path, _ in self.rc_graph_merger_list:
        #         self.send_log(rc, script_path)

        elapsed_time = end_time - start_time
        err_print(f"Execution took {elapsed_time} seconds")

        if not self.kill_target and hasattr(self, "first_request_time"):
            elapsed_time_for_killing = end_time - self.first_request_time
            self.worker.last_exec_time_dict[self.script_name] = elapsed_time_for_killing * 1000
            err_print(f"Updating last execution time for \"{self.script_name}\" to {elapsed_time_for_killing} seconds, target is \"{self.kill_target}\"")
        else:
            err_print(f"Not updating last execution time, target is \"{self.kill_target}\"")

        if self.debug:
            result1 = subprocess.run(['du', '-h', '-d0', config.PASH_TMP_PREFIX], capture_output=True, text=True, check=True)
            if self.ft == "optimized":
                result2 = subprocess.run(['du', '-h', '-d0', self.fish_out_prefix], capture_output=True, text=True, check=True)
                # Fish outs are inside tempdirs
                err_print(f"Temp dir size | Fish out size: {result1.stdout.split()[0]} | {result2.stdout.split()[0]}")
            else:
                err_print(f"Temp dir size: {result1.stdout.split()[0]}")

        # Remove directory (except for unix socket files)
        # with os.scandir(self.pash_tmp_prefix) as entries:
        #     for entry in entries:
        #         mode = os.stat(entry.path).st_mode
        #         isSocket = stat.S_ISSOCK(mode)
        #         err_print(entry.path, isSocket)
        #         if isSocket:
        #             err_print(f"Skipping unix socket file {entry.path}")
        #             continue
        #         else:
        #             if entry.is_dir():
        #                 shutil.rmtree(entry)
        #             else:
        #                 if entry.name == "dspash_socket":
        #                     continue
        #                 else:
        #                     os.remove(entry.path)
        
        shutil.rmtree(self.pash_tmp_prefix)
        err_print(f"Temporary directory deleted: {self.pash_tmp_prefix}")

    def send_log(self, rc: subprocess.Popen, script_path: str):
        name = self.debug['name']
        url = self.debug['url']
        with open(script_path, 'r') as file:
            shell_script = file.read()

        try:
            # timeout is set to 10s for debuggin
            _, err = rc.communicate(timeout=10)
        except:
            err_print("process timedout")
            rc.kill()
            _, err = rc.communicate()

        if err is None:
            err = b""

        response = {
            'name': name,
            'returncode': rc.returncode,
            'stderr': err.decode("UTF-8"),
            'shellscript': shell_script,
        }

        requests.post(url=url, json=response)

    def handle_setup_request(self):
        self.debug = int(self.request['debug'])
        self.pool_size = int(self.request['pool_size'])
        self.ft = self.request['ft']
        self.script_name = self.request['script_name']
        # here we won't be receiving the address but we are still required
        # to store it to not update last exec times
        self.kill_target = self.request.get('kill_target', None)

        global DEBUG
        if self.debug:
            DEBUG = True
            err_print(f"Debugging enabled, level: {self.debug}")
        else:
            err_print(f"Debugging disabled, level: {self.debug}")
            DEBUG = False

        if self.debug > 1:
            self.stderr=None
        else:
            self.stderr=subprocess.DEVNULL

        if self.ft == "optimized":
            self.event_loop = EventLoop(self)
            err_print(f"Starting {self.event_loop.name} with pool size {self.pool_size}")
            self.event_loop.start()
        
        # Run mktemp -d to create a unique temporary directory
        result = subprocess.run(['mktemp', '-d', '/tmp/pash_XXXXXXX'], capture_output=True, text=True)

        # Extract the directory path from the output
        if result.returncode == 0:
            temp_directory = result.stdout.strip()
            # Set this path as an environment variable
            os.environ['PASH_TMP_PREFIX'] = temp_directory
            config.PASH_TMP_PREFIX = temp_directory
            self.pash_tmp_prefix = temp_directory
            err_print(f"Temporary directory created: {temp_directory}")
        else:
            raise Exception(f"Failed to create temporary directory: {result.stderr}")
        
        datastream_directory = os.path.join(temp_directory, 'datastream')
        os.makedirs(datastream_directory, exist_ok=True)
        os.environ['FISH_OUT_PREFIX'] = datastream_directory
        self.env_var['FISH_OUT_PREFIX'] = datastream_directory
        self.fish_out_prefix = datastream_directory
        err_print(f"Datastream directory created: {datastream_directory}")

    def handle_kill_node(self):
        self.kill_target = self.request['kill_target']
        if self.kill_target == HOST:
            if self.request['kill_delay']:
                delay = int(self.request['kill_delay'])
            else:
                delay = self.worker.last_exec_time_dict[self.script_name] / 2
            err_print(f"Starting killer thread for script \"{self.script_name}\"! Delay is {delay} millis.")
            Thread(target=kill, args=(delay, self.event_loop)).start()
        else:
            # This should never happen with the new approach
            err_print(f"Kill target {self.kill_target} is not this node, skipping kill")

    def handle_exec_graph_request(self):
        if self.ft != "optimized":
            save_configs(self.request['graph'], self.dfs_configs_paths)

        if self.first_request:
            self.first_request_time = time.time()
            self.first_request = False

        config.config['shell_variables'] = self.request['shell_variables']
        if self.debug:
            add_debug_flags(self.request['graph'])
            stderr = subprocess.PIPE
        else:
            stderr = None

        script_path = to_shell_file(self.request['graph'], config.pash_args)

        e = os.environ.copy()
        e['PASH_TOP'] = PASH_TOP
        e['DISH_TOP'] = DISH_TOP

        # store functions
        functions_file = create_filename(dir=self.pash_tmp_prefix, prefix='pashFuncs')
        write_file(functions_file, self.request['functions'])
        cmd = f"source {functions_file}; source {script_path}"
        err_print(f"executing {cmd}")

        rc = subprocess.Popen(cmd, env=e, executable="/bin/bash", shell=True, stderr=self.stderr, preexec_fn=os.setsid)
        self.rc_graph_merger_list.append((rc, script_path, self.request['merger_id']))

    def handle_batch_exec_graph(self):
        # this in None for now anyways
        # config.config['shell_variables'] = self.request['shell_variables']

        err_print(f"Batch exec graph request: {len(self.request['regulars'])} regulars and {len(self.request['mergers'])} mergers")
        err_print(f"PASH_TMP_PREFIX {self.pash_tmp_prefix}")

        if self.first_request:
            self.first_request_time = time.time()
            self.first_request = False

        functions_file = create_filename(dir=self.pash_tmp_prefix, prefix='pashFuncs')
        write_file(functions_file, self.request['functions'])

        for regular in self.request['regulars']:
            if self.debug > 1:
                add_debug_flags(regular)
            script_path = to_shell_file(regular, config.pash_args)
            cmd = f"source {functions_file}; source {script_path}"
            self.event_loop.regular_queue.append((cmd, script_path, self.request['merger_id']))

        for merger in self.request['mergers']:
            if self.debug > 1:
                add_debug_flags(merger)
            script_path = to_shell_file(merger, config.pash_args)
            cmd = f"source {functions_file}; source {script_path}"
            if self.event_loop.is_alive():
                err_print(f"Executing merger {cmd}")
                rc = subprocess.Popen(cmd, env=self.env_var, executable="/bin/bash", shell=True, stderr=self.stderr, preexec_fn=os.setsid)
                self.rc_graph_merger_list.append((rc, script_path, self.request['merger_id']))

    def handle_kill_subgraphs_request(self):
        merger_id = self.request['merger_id']
        # Record the start time
        start_time = time.time()
        err_print(f"Killing subgraphs, merger_id: {merger_id}")

        # Try to terminate the subprocess gracefully
        kill_count = 0
        for rc, _, m in self.rc_graph_merger_list:
            if merger_id != -1 and m != merger_id:
                continue
            os.killpg(os.getpgid(rc.pid), signal.SIGKILL)
            kill_count += 1

        for rc, _, m in self.rc_graph_merger_list:
            if merger_id != -1 and m != merger_id:
                continue

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

        err_print(f"Killed subgraphs in {elapsed_time} seconds, {kill_count} subgraphs killed")
        return {}

    def handle_resurrect_request(self):
        global DEBUG
        DEBUG_TMP = DEBUG
        # Temporarily enable DEBUG state
        DEBUG = True
        script_path = "$DISH_TOP/docker-hadoop/datanode/run.sh"
        subprocess.run("/bin/bash " + script_path + " --resurrect", shell=True)
        time.sleep(1)
        err_print("Just brought the current node back up!")
        # Restore DEBUG state
        DEBUG = DEBUG_TMP


class EventLoop(Thread):
    def __init__(self, handler: RequestHandler):
        Thread.__init__(self)
        self.name = "Eventloop-" + handler.name
        self.handler = handler
        self.quit = Event()
        self.regular_queue = []
        self.running_processes: List[subprocess.Popen] = []

    def run(self):
        while not self.quit.is_set():
            # Remove finished processes
            for i in range(len(self.running_processes)-1, -1, -1):
                if self.running_processes[i].poll() is not None:
                    self.running_processes.pop(i)

            # Run until pool is full
            while self.regular_queue and len(self.running_processes) < self.handler.pool_size:
                cmd, script_path, merger_id = self.regular_queue.pop()
                # err_print(f"executing {cmd}")
                rc = subprocess.Popen(cmd, env=self.handler.env_var, executable="/bin/bash", shell=True, stderr=self.handler.stderr, preexec_fn=os.setsid)
                # rc = subprocess.Popen(f"time $({cmd})", env=self.handler.env_var, executable="/bin/bash", shell=True, stderr=None, preexec_fn=os.setsid)
                self.running_processes.append(rc)
                self.handler.rc_graph_merger_list.append((rc, script_path, merger_id))

            # Sleep for a bit to not busy wait
            self.quit.wait(0.1)


def err_print(*args):
    if DEBUG:
        print(*args, file=sys.stderr, flush=True)


def send_success(conn, body, msg=""):
    request = {
        'status': 'OK',
        'body': body,
        'msg': msg
    }
    send_msg(conn, encode_request(request))


def kill(delay: int, event_loop: EventLoop):
    err_print(f"Will kill after delay for {delay}")
    t0 = time.time()
    time.sleep(delay / 1000)
    t1 = time.time()
    err_print(f"Killing now {t1 - t0} seconds after delay")
    if event_loop:
        event_loop.quit.set()
        event_loop.join()
        t2 = time.time()
        err_print(f"Event loop joined {t2 - t1} seconds after")

    script_path = "$DISH_TOP/runtime/scripts/killall.sh"
    subprocess.run("/bin/sh " + script_path, shell=True)
    time.sleep(1)
    subprocess.run("/bin/sh " + script_path, shell=True)
    t3 = time.time()
    err_print(f"Just killed! {t3 - t2} seconds after")


def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    config.add_common_arguments(parser)
    # Add an additional "resurrect" arg for bringing node back up. 
    # This is not added to config.py because this is an internal flag
    parser.add_argument('--resurrect', action='store_true', help='To bring up the current node')
    args = parser.parse_args()
    config.pash_args = args
    # Initialize the log file
    config.init_log_file()
    if not config.config:
        config.load_config(args.config_path)


def init():
    parse_args()
    config.LOGGING_PREFIX = f"Worker {PORT}: "
    config.annotations = load_annotation_files(config.config['distr_planner']['annotations_dir'])
    pash_runtime.runtime_config = config.config['distr_planner']
    # needed because graphs could have multiples sinks
    config.pash_args.termination = ""


def main():
    init()
    worker = Worker()
    worker.run()


if __name__ == "__main__":
    main()
