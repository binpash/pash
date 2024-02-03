import socket
import os
import time
import queue
import pickle
import json
from typing import List

from dspash.socket_utils import SocketManager, encode_request, decode_request, send_msg, recv_msg
from definitions.ir.file_id import FileId
from util import log
from dspash.ir_helper import prepare_graph_for_remote_exec, to_shell_file
from dspash.utils import read_file
import config 
import copy
import requests

PORT = 65425        # Port to listen on (non-privileged ports are > 1023)
# TODO: get the url from the environment or config
DEBUG_URL = f'http://{socket.getfqdn()}:5001' 

class WorkerConnection:
    def __init__(self, name, host, port):
        self.name = name
        self._host = socket.gethostbyaddr(host)[2][0] # get ip address in case host needs resolving
        self._port = port
        self._running_processes = 0
        self._online = True
        # assume client service is running, can add a way to activate later
        try:
            self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._socket.connect((self._host, self._port))
        except Exception as e:
            log(f"Failed to connect to {self._host}:{self._port} with error {e}")
            self._online = False
        
    def is_online(self):
        # TODO: create a ping to confirm is online
        return self._online

    def get_running_processes(self):
        # request_dict = { 'type': 'query',
        #                 'fields': ['running_processes']
        # }
        # with self.socket.connect((self.host, self.port)):
        #     self.socket.send(request)
        #     # TODO wait until the command exec finishes and run this in parallel?
        #     answer = self.socket.recv(1024)
        return self._running_processes

    def send_graph_exec_request(self, graph, shell_vars, functions, args, worker_timeout=0) -> bool:
        request_dict = { 'type': 'Exec-Graph',
                        'graph': graph,
                        'functions': functions,
                        'shell_variables': None, # Doesn't seem needed for now
                        'debug': None,
                        'worker_timeout': worker_timeout,
                        'kill': args.kill,
                    }
        if args.debug:
            request_dict['debug'] = {'name': self.name, 'url': f'{DEBUG_URL}/putlog'}

        request = encode_request(request_dict)
        #TODO: do I need to open and close connection?
        send_msg(self._socket, request)
        # TODO wait until the command exec finishes and run this in parallel?
        retries = 0
        MAX_RETRIES = 2
        RETRY_DELAY = 1
        self._socket.settimeout(5)
        response_data = None
        response = None
        while retries < MAX_RETRIES and not response_data:
            try:
                response_data = recv_msg(self._socket)
                response = decode_request(response_data)
                log(f"Socket {self._socket.getsockname()} recieved response from {self._socket.getpeername()}, response: {response}")
            except socket.timeout:
                log(f"Timeout encountered. Retry {retries + 1} of {MAX_RETRIES}.")
                retries += 1
                time.sleep(RETRY_DELAY)
        if not response_data or response['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
        else:
            # self._running_processes += 1 #TODO: decrease in case of failure or process ended
            return True

    def close(self):
        self._socket.send(encode_request({"type": "Done"}))
        self._socket.close()
        self._online = False

    def _wait_ack(self):
        confirmation = self._socket.recv(4096)
        if not confirmation or decode_request(confirmation).status != "OK":
            return False
            # log(f"Confirmation not recieved {confirmation}")
        else:
            return True

    def __str__(self):
        return f"Worker {self._host}:{self._port}"

    def host(self):
        return self._host

class WorkersManager():
    def __init__(self, workers: List[WorkerConnection] = []):
        self.workers = workers
        self.host = socket.gethostbyname(socket.gethostname())
        self.args = copy.copy(config.pash_args)
        # Required to create a correct multi sink graph
        self.args.termination = "" 

    def get_worker(self, fids: List[FileId] = None) -> WorkerConnection:
        if not fids:
            fids = []

        best_worker = None  # Online worker with least work
        for worker in self.workers:
            if not worker.is_online():
                log(f"Skipping {worker} because it is offline")
                continue
            
            # Skip if any provided fid isn't available on the worker machine
            if any(map(lambda fid: not fid.is_available_on(worker.host()), fids)):
                continue

            if best_worker is None or best_worker.get_running_processes() > worker.get_running_processes():
                best_worker = worker

        if best_worker == None:
            raise Exception("no workers online where the data is stored")

        return best_worker

    def add_worker(self, name, host, port):
        self.workers.append(WorkerConnection(name, host, port))

    def add_workers_from_cluster_config(self, config_path):
        with open(config_path, 'r') as f:
            cluster_config = json.load(f)

        workers = cluster_config["workers"]
        for name, worker in workers.items():
            host = worker['host']
            port = worker['port']
            self.add_worker(name, host, port)
            
            
    def run(self):
        workers_manager = self
        workers_manager.add_workers_from_cluster_config(os.path.join(config.PASH_TOP, 'cluster.json'))
        if workers_manager.args.debug:
            try:
                requests.post(f'{DEBUG_URL}/clearall') # clears all the debug server logs
            except Exception as e:
                log(f"Failed to connect to debug server with error {e}\n")
                workers_manager.args.debug = False # Turn off debugging

        dspash_socket = SocketManager(os.getenv('DSPASH_SOCKET'))
        while True:
            request, conn = dspash_socket.get_next_cmd()
            if request.startswith("Done"):
                dspash_socket.close()
                break
            elif request.startswith("Exec-Graph"):
                args = request.split(':', 1)[1].strip()
                filename, declared_functions_file = args.split()
                numExecutedSubgraphs = 0
                numTotalSubgraphs = None
                crashed_worker = workers_manager.args.worker_timeout_choice if workers_manager.args.worker_timeout_choice != '' else "worker1" # default to be worker1
                try:
                    while not numTotalSubgraphs or numExecutedSubgraphs < numTotalSubgraphs:
                        # In the naive fault tolerance, we want all workers to receive its subgraph(s) without crashing
                        # if a crash happens, we'll re-split the IR and do it again until scheduling is done without any crash.
                        numExecutedSubgraphs = 0
                        worker_subgraph_pairs, shell_vars, main_graph = prepare_graph_for_remote_exec(filename, workers_manager.get_worker)
                        if numTotalSubgraphs == None:
                            numTotalSubgraphs = len(worker_subgraph_pairs)
                        script_fname = to_shell_file(main_graph, workers_manager.args)
                        log("Master node graph stored in ", script_fname)

                        # Read functions
                        log("Functions stored in ", declared_functions_file)
                        declared_functions = read_file(declared_functions_file)

                        # Execute subgraphs on workers
                        for worker, subgraph in worker_subgraph_pairs:
                            worker: WorkerConnection
                            worker_timeout = workers_manager.args.worker_timeout if worker.name == crashed_worker and workers_manager.args.worker_timeout else 0
                            
                            try:
                                worker.send_graph_exec_request(subgraph, shell_vars, declared_functions, workers_manager.args, worker_timeout)
                                numExecutedSubgraphs += 1
                            except Exception as e:
                                # worker timeout
                                worker.close()
                                log(f"{worker} closed")
                    # Report to main shell a script to execute
                    # Delay this to the very end when every worker has received the subgraph
                    response_msg = f"OK {script_fname}"
                    dspash_socket.respond(response_msg, conn)
                except Exception as e:
                    print(e)

                    

            else:
                raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    WorkersManager().run()
