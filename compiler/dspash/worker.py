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
from threading import Event, Thread

DISH_TOP = os.environ['DISH_TOP'] 
PASH_TOP = os.environ['PASH_TOP']
sys.path.append(os.path.join(PASH_TOP, "compiler"))

from dspash.utils import create_filename, write_file
from dspash.ir_helper import save_configs, to_shell_file, add_debug_flags
from dspash.socket_utils import send_msg, recv_msg
import pash_runtime
from annotations import load_annotation_files
import config

import cProfile
import pstats
import io

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
        self.killed = False
        self.script_name = ""

    def run(self):
        with self.conn:
            err_print('Connected by', self.addr, 'name', self.name)
            self.start_time = time.time()
            doCleanup = False

            while True:
                data = recv_msg(self.conn)
                if not data:
                    break
                self.request = decode_request(data)
                self.rh_print("Got a new request", self.request['type'])

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
                    self.rh_print(f"Unsupported request {self.request}")
                send_success(self.conn, self.body)
            if doCleanup:
                self.cleanup()
        if self.debug > 2:
            self.profiler.disable()
            s = io.StringIO()
            ps = pstats.Stats(self.profiler, stream=s).sort_stats('cumulative')
            ps.print_stats()
            with open(f"W_profile.log", "w") as f:
                f.write(s.getvalue())

    def handle_setup_request(self):
        self.debug = int(self.request['debug'])
        if self.request['pool_size'] == "auto":
            self.pool_size = os.cpu_count()
        else:
            self.pool_size = int(self.request['pool_size'])
        self.ft = self.request['ft']
        self.script_name = self.request['script_name']
        # here we won't be receiving the address but we are still required
        # to store it to not update last exec times
        self.kill_target = self.request.get('kill_target', None)

        global DEBUG
        if self.debug:
            DEBUG = True
            self.rh_print(f"Debugging enabled, level: {self.debug}")
        else:
            self.rh_print(f"Debugging disabled, level: {self.debug}")
            DEBUG = False

        if self.debug > 1:
            self.stderr=None
        else:
            self.stderr=subprocess.DEVNULL

        if self.debug > 2:
            self.profiler = cProfile.Profile()
            self.profiler.enable()

        if self.ft == "optimized" or self.ft == "dynamic":
            self.event_loop = EventLoop(self)
            self.rh_print(f"Starting {self.event_loop.name} with pool size {self.pool_size}")
            self.event_loop.start()
        
        self.time_recorder = TimeRecorder(self)
        self.rh_print(f"Starting {self.time_recorder.name}")
        self.time_recorder.start()
        
        # Run mktemp -d to create a unique temporary directory
        result = subprocess.run(['mktemp', '-d', '/tmp/pash_XXXXXXX'], capture_output=True, text=True)

        # Extract the directory path from the output
        if result.returncode == 0:
            temp_directory = result.stdout.strip()
            # Set this path as an environment variable
            os.environ['PASH_TMP_PREFIX'] = temp_directory
            config.PASH_TMP_PREFIX = temp_directory
            self.pash_tmp_prefix = temp_directory
            self.rh_print(f"Temporary directory created: {temp_directory}")
        else:
            raise Exception(f"Failed to create temporary directory: {result.stderr}")
        
        datastream_directory = os.path.join(temp_directory, 'datastream')
        os.makedirs(datastream_directory, exist_ok=True)
        os.environ['FISH_OUT_PREFIX'] = datastream_directory
        self.env_var['FISH_OUT_PREFIX'] = datastream_directory
        self.fish_out_prefix = datastream_directory
        self.rh_print(f"Datastream directory created: {datastream_directory}")

    def handle_kill_node(self):
        self.kill_target = self.request['kill_target']
        if self.kill_target == HOST:
            if self.request['kill_delay']:
                # We take delay input in millis, so convert to seconds
                delay = float(self.request['kill_delay']) / 1000
            else:
                if self.script_name in self.worker.last_exec_time_dict:
                    delay = self.worker.last_exec_time_dict[self.script_name] / 2
                else:
                    # This can happen if no exec request was received due to small input or, if you forgot to run faultless version first.
                    # It makes sense to kill here with a small delay to not cause problems for resurrect
                    delay = 1
                    self.rh_print(f"Script don't have any last_exec_time")
            self.rh_print(f"Starting killer thread, delay is {delay}")
            Thread(target=self.kill, args=(delay,)).start()
        else:
            # This should never happen with the new approach
            self.rh_print(f"Kill target {self.kill_target} is not this node, skipping kill")

    def kill(self, delay: int):
        self.rh_print("KT: Started")
        time.sleep(delay)
        self.rh_print("KT: Killing now")
        self.killed = True

        if self.event_loop:
            self.event_loop.quit.set()
            self.event_loop.join()
            self.rh_print("KT: Event loop joined")

        script_path = "$DISH_TOP/runtime/scripts/killall.sh"
        subprocess.run("/bin/sh " + script_path, shell=True)
        time.sleep(1)

        self.rh_print("KT: Killing once more")
        subprocess.run("/bin/sh " + script_path, shell=True)
        self.rh_print("KT: Kill finished")

    def handle_exec_graph_request(self):
        self.rh_print(f"Exec graph request received, id: {self.request['graph'].id}, merger {self.request['graph'].merger}")
        if self.killed:
            self.rh_print(f"Killed, skipping exec request, id: {self.request['graph'].id}, merger {self.request['graph'].merger}")
            return

        if self.ft != "optimized":
            save_configs(self.request['graph'], self.dfs_configs_paths)

        config.config['shell_variables'] = self.request['shell_variables']
        if self.debug > 1:
            add_debug_flags(self.request['graph'])

        script_path = to_shell_file(self.request['graph'], config.pash_args)

        e = os.environ.copy()
        e['PASH_TOP'] = PASH_TOP
        e['DISH_TOP'] = DISH_TOP

        # store functions
        functions_file = create_filename(dir=self.pash_tmp_prefix, prefix='pashFuncs')
        write_file(functions_file, self.request['functions'])
        cmd = f"source {functions_file}; source {script_path}"
        # self.rh_print(f"executing {cmd}")

        rc = subprocess.Popen(cmd, env=e, executable="/bin/bash", shell=True, stderr=self.stderr, preexec_fn=os.setsid)
        self.rc_graph_merger_list.append((rc, script_path, self.request['merger_id']))

        # self.rh_print(f"Exec graph request handled, id: {self.request['graph'].id}, merger {self.request['graph'].merger}")

    def handle_batch_exec_graph(self):
        if self.killed:
            self.rh_print("Killed, skipping exec request")
            return
        # this in None for now anyways
        # config.config['shell_variables'] = self.request['shell_variables']

        self.rh_print(f"Batch exec graph request: {len(self.request['regulars'])} regulars and {len(self.request['mergers'])} mergers")
        self.rh_print(f"PASH_TMP_PREFIX {self.pash_tmp_prefix}")

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
                self.rh_print(f"Executing merger {cmd}")
                rc = subprocess.Popen(cmd, env=self.env_var, executable="/bin/bash", shell=True, stderr=self.stderr, preexec_fn=os.setsid)
                self.rc_graph_merger_list.append((rc, script_path, self.request['merger_id']))

    def __kill_pg(self, rc):
            try:
                pg_id = os.getpgid(rc.pid)  # Get the process group ID
                os.killpg(pg_id, signal.SIGKILL)
                self.rh_print(f"Process group {pg_id} killed successfully.")
            except ProcessLookupError:
                self.rh_print(f"No such process group for PID {rc.pid}. It may have already terminated.")

    def handle_kill_subgraphs_request(self):
        merger_id = self.request['merger_id']
        # Record the start time
        start_time = time.time()
        self.rh_print(f"Killing subgraphs, merger_id: {merger_id}")

        # If eventloop is enabled, mark the merger for removal from the eventloop queue
        if self.event_loop:
            self.event_loop.merger_to_remove = merger_id
            self.rh_print(f"Merger {merger_id} marked for removal for eventloop queue")
            while self.event_loop.merger_to_remove != -2:
                self.rh_print(f"Waiting for merger {merger_id} to be removed from eventloop queue")
                self.event_loop.quit.wait(0.1)
            self.rh_print(f"Merger {merger_id} removed from eventloop queue")

        # Try to terminate the subprocess gracefully
        kill_count = 0
        for rc, _, m in self.rc_graph_merger_list:
            if merger_id != -1 and m != merger_id:
                continue
            self.__kill_pg(rc)
            kill_count += 1

        for rc, _, m in self.rc_graph_merger_list:
            if merger_id != -1 and m != merger_id:
                continue

            try:
                # Wait for the subprocess to terminate, with a timeout
                rc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                # If the subprocess is still running after the timeout, kill it
                self.rh_print(f"Terminating subprocess {rc.pid} failed, killing it forcefully")
                self.__kill_pg(rc)
                rc.wait()

        # Record the end time
        end_time = time.time()

        # Calculate the time elapsed
        elapsed_time = end_time - start_time

        self.rh_print(f"Killed subgraphs in {elapsed_time} seconds, {kill_count} subgraphs killed")
        return {}

    def handle_resurrect_request(self):
        global DEBUG
        DEBUG_TMP = DEBUG
        # Temporarily enable DEBUG state
        DEBUG = True
        script_path = "$DISH_TOP/docker-hadoop/datanode/run.sh"
        subprocess.run("/bin/bash " + script_path + " --resurrect", shell=True)
        time.sleep(1)
        self.rh_print("Just brought the current node back up!")
        # Restore DEBUG state
        DEBUG = DEBUG_TMP

    def cleanup(self):
        self.rh_print("Connection ended. Subprocess count:", len(self.rc_graph_merger_list))

        for rc, _, _ in self.rc_graph_merger_list:
            rc.wait()

        if self.ft == "optimized" or self.ft == "dynamic":
            self.event_loop.quit.set()
            self.event_loop.join()
            self.rh_print(f"{self.event_loop.name} joined")

        self.time_recorder.quit.set()
        self.time_recorder.join()
        end_time = self.time_recorder.end_time
        self.rh_print(f"{self.time_recorder.name} joined")

        elapsed_time = end_time - self.start_time
        self.rh_print(f"Execution took {elapsed_time} seconds for \"{self.script_name}\"")

        if not self.kill_target:
            self.worker.last_exec_time_dict[self.script_name] = elapsed_time
            self.rh_print(f"Updating last execution time to {elapsed_time} seconds")
        else:
            self.rh_print(f"Not updating last execution time")

        if self.debug:
            result1 = subprocess.run(['du', '-h', '-d0', config.PASH_TMP_PREFIX], capture_output=True, text=True, check=True)
            result2 = subprocess.run(['du', '-h', '-d0', self.fish_out_prefix], capture_output=True, text=True, check=True)
            # Fish outs are inside tempdirs
            self.rh_print(f"Temp dir size | Fish out size: {result1.stdout.split()[0]} | {result2.stdout.split()[0]}")

        if self.debug > 1:
            self.rh_print(f"Skipping deleting temporary directory ({self.pash_tmp_prefix}) as debug enabled")
            datastream_directory = os.path.join(self.pash_tmp_prefix, 'datastream')
            if os.path.exists(datastream_directory):
                shutil.rmtree(datastream_directory)
                self.rh_print(f"Datastream directory deleted: {datastream_directory}")
        else:
            if os.path.exists(self.pash_tmp_prefix):
                shutil.rmtree(self.pash_tmp_prefix)
                self.rh_print(f"Temporary directory deleted: {self.pash_tmp_prefix}")

        self.rh_print("Cleanup done")
        self.rh_print("=" * 100)
        self.rh_print("=" * 100)
        self.rh_print("=" * 100)
        self.rh_print("=" * 100)
        self.rh_print("=" * 100)

    def rh_print(self, *args):
        if DEBUG:
            err_print(f"{self.name} \"{self.ft}\" \"{self.kill_target}\" \"{self.script_name}\" {(time.time() - self.start_time):10.6f}", *args)


class EventLoop(Thread):
    def __init__(self, handler: RequestHandler):
        Thread.__init__(self)
        self.name = "EventLoop-" + handler.name
        self.handler = handler
        self.quit = Event()
        self.regular_queue = []
        self.running_processes: List[subprocess.Popen] = []
        # -2 means nothing to remove, -1 remove all, else remove merger_id and it's dependencies
        self.merger_to_remove = -2

    def run(self):
        while not self.quit.is_set():
            # Remove finished processes
            for i in range(len(self.running_processes)-1, -1, -1):
                if self.running_processes[i].poll() is not None:
                    self.running_processes.pop(i)

            self.handler.rh_print(f"Eventloop: Current length {len(self.regular_queue)}, running {len(self.running_processes)}, to_remove {self.merger_to_remove}")
            if self.merger_to_remove != -2:
                # It's guaranteed by orchestrator that we will never recevie a new request until we finish this operation
                self.handler.rh_print(f"Eventloop: Merger {self.merger_to_remove} received, current length {len(self.regular_queue)}")
                if self.merger_to_remove == -1:
                    self.regular_queue.clear()
                else:
                    # t = (cmd, script_path, merger_id)
                    self.regular_queue = [t for t in self.regular_queue if t[2] != self.merger_to_remove]
                self.merger_to_remove = -2
                self.handler.rh_print(f"Eventloop: Merger {self.merger_to_remove} removed, current length {len(self.regular_queue)}")

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


class TimeRecorder(Thread):
    def __init__(self, handler: RequestHandler):
        Thread.__init__(self)
        self.name = "TimeRecorder-" + handler.name
        self.handler = handler
        self.quit = Event()
        self.last_len = 0
        self.end_time = 0

    def run(self):
        while not self.quit.is_set():
            curr_len = len(self.handler.rc_graph_merger_list)
            if curr_len != self.last_len:
                for rc, _, _ in self.handler.rc_graph_merger_list[:]:
                    rc.wait()
                self.last_len = curr_len
                self.end_time = time.time()
                err_print(f"{(self.end_time - self.handler.start_time):10.6f} TimeRecorder: end_time recorded {self.end_time}, current length {curr_len}")

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
