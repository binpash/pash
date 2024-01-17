import socket
import os
import time
import queue
import pickle
import json
import collections

from dspash.socket_utils import SocketManager, encode_request, decode_request, send_msg, recv_msg
from util import log
from dspash.ir_helper import prepare_graph_for_remote_exec, to_shell_file, get_best_worker_for_subgraph, get_remote_pipe_update_candidates, update_subgraphs, get_worker_subgraph_map, update_remote_pipe_addr
from dspash.utils import read_file
import config 
import copy
import requests
import definitions.ir.nodes.remote_pipe as remote_pipe
from ir_to_ast import to_shell


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

    def send_graph_exec_request(self, graph, shell_vars, functions, debug=False, worker_timeout=0) -> bool:
        request_dict = { 'type': 'Exec-Graph',
                        'graph': graph,
                        'functions': functions,
                        'shell_variables': None, # Doesn't seem needed for now
                        'debug': None,
                        'worker_timeout': worker_timeout    
                    }
        if debug:
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
        
    def abortAll(self, debug=False):
        request_dict = { 'type': 'abortAll',
                        'debug': None
                    }
        if debug:
            request_dict['debug'] = {'name': self.name, 'url': f'{DEBUG_URL}/debug'}
        request = encode_request(request_dict)
        send_msg(self._socket, request)
        response_data = recv_msg(self._socket)
        if not response_data or decode_request(response_data)['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
        else:
            response = decode_request(response_data)
            print(response)

    def getDiscoveryServerLog(self, debug=False):
        request_dict = { 'type': 'getDiscoveryServerLog',
                        'debug': None
                    }
        if debug:
            request_dict['debug'] = {'name': self.name, 'url': f'{DEBUG_URL}/debug'}

        request = encode_request(request_dict)
        send_msg(self._socket, request)

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
    def __init__(self, workers: WorkerConnection = []):
        self.workers = workers
        self.host = socket.gethostbyname(socket.gethostname())
        self.args = copy.copy(config.pash_args)
        # Required to create a correct multi sink graph
        self.args.termination = "" 

    def get_worker(self, fids = None) -> WorkerConnection:
        if not fids:
            fids = []

        best_worker = None  # Online worker with least work
        for worker in self.workers:
            if not worker.is_online():
                continue
            
            # Skip if any provided fid isn't available on the worker machine
            if any(map(lambda fid: not fid.is_available_on(worker.host()), fids)):
                continue

            if best_worker is None or best_worker.get_running_processes() > worker.get_running_processes():
                best_worker = worker

        if best_worker == None:
            raise Exception("no workers online where the date is stored")

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
            

    def deep_copy_worker_subgraph_pairs(self, worker_subgraph_pairs):
        # worker is not serializable and shallow copy is fine
        deep_copy = []
        for worker, subgraph in worker_subgraph_pairs:
            deep_copy.append([worker, copy.deepcopy(subgraph)])
        return collections.deque(deep_copy)

            
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

                crashed_worker_choice = workers_manager.args.worker_timeout_choice if workers_manager.args.worker_timeout_choice != '' else "worker1" # default to be worker1
                try:
                    # In the naive fault tolerance, we want all workers to receive its subgraph(s) without crashing
                    # if a crash happens, we'll re-split the IR and do it again until scheduling is done without any crash.
                    worker_subgraph_pairs, shell_vars, main_graph, file_id_gen, input_fifo_map = prepare_graph_for_remote_exec(filename, workers_manager.get_worker)
                    
                    # Read functions
                    log("Functions stored in ", declared_functions_file)
                    declared_functions = read_file(declared_functions_file)

                    # Execute subgraphs on workers
                    # worker_subgraph_map = get_worker_subgraph_map(worker_subgraph_pairs)
                    worker_subgraph_pairs = collections.deque(worker_subgraph_pairs)
                    #TODO: deepcopy
                    worker_subgraph_pairs_backup = self.deep_copy_worker_subgraph_pairs(worker_subgraph_pairs)
                    

                    while worker_subgraph_pairs:
                        worker, subgraph = worker_subgraph_pairs.popleft()
                        print(worker.name, worker)
                        print(to_shell(subgraph, workers_manager.args))
                        print('---------------------------------')
                        worker_timeout = workers_manager.args.worker_timeout if worker.name == crashed_worker_choice and workers_manager.args.worker_timeout else 0
                        try:
                            worker.send_graph_exec_request(subgraph, shell_vars, declared_functions, workers_manager.args.debug, worker_timeout)
                        except Exception as e:
                            # For each subgraph assigned to worker (now crashed), do:
                            #   1) find the next best worker as the replacement_worker
                            #   2) get every remote_pipe (call it RP) in the subgraph whose addr is the crashed_worker's host
                            #       and all neighboring remote_pipes (call it RP-N) to RP
                            #       which are remote_pipes on other subgraphs that communicated with RP.
                            #       we will need to update addresses for RP | RP-N from the crashed_worker's host to replacement_worker's host
                            #       Note: we don't want to update host for evert remote_pipe on the subgraph because it could have been
                            #              a remote_read pipe reading from a healthy worker. In this case we want to keep it as it is
                            #   3) Update addr for all remote_write pipes on RP and all remote_read pipes on RP-N from the crashed_worker's host to replacement_worker's host
                            #       (we want to hold off updating addr for remote_read pipe on RP as subgraph1 and subgraphs2 can be assigned to different workers)                                    
                            #   *Note: subgraphs assigned to the same crashed worker may be assigned to different replacement workers
                            #   For example, 
                            #   worker1 at addr 5 has subgraph1 and subgraph2, and subgraph1 is a merger node such that
                            #           subgraph1 has a remote_read 5 uuid123 and also a remote_write 5 uuid567
                            #           subgraph2 has a remote_write 5 uuid123
                            #   When worker1 crashes, subgraph1 can be assigned to worker2 at addr 6 and subgraph2 to worker3
                            #
                            #
            
                            # Shuts down crashed worker gracefully
                            crashed_worker = worker
                            if crashed_worker.is_online():
                                crashed_worker.close()
                                log(f"{crashed_worker} closed")
                            else:
                                log(f"{crashed_worker} is already closed")

                            # TODO: do we want to send all relevant workers a msg to abort current subgraphs?
                            # This is buggy so not used
                            for worker in self.workers:
                                if worker.is_online():
                                    worker.abortAll(workers_manager.args.debug)

                            subgraphs = [pair[1] for pair in worker_subgraph_pairs_backup]
                            for i in range(len(worker_subgraph_pairs_backup)):
                                worker, subgraph = worker_subgraph_pairs_backup[i]
                                if worker == crashed_worker:
                                    # Step 1
                                    # find the next best worker for the subgraph and push (worker, subgraph) pair back to the deque
                                    replacement_worker = get_best_worker_for_subgraph(subgraph, workers_manager.get_worker)
                                    worker_subgraph_pairs_backup[i][0] = replacement_worker

                                    # Step 2
                                    # Find all write RP on subgraph and all read RP-N on neighboring subgraphs (including main_graph)
                                    update_candidates = get_remote_pipe_update_candidates(subgraphs + [main_graph], subgraph, crashed_worker)
                                    for candidate in update_candidates:
                                        print(candidate.get_uuid(), candidate.get_host())
                                    # Step 3
                                    for update_candidate in update_candidates:
                                        update_remote_pipe_addr(update_candidate, replacement_worker.host())
                                
                            
                            # Update meta-data
                                        # TODO: DEEPCOPY
                            worker_subgraph_pairs = self.deep_copy_worker_subgraph_pairs(worker_subgraph_pairs_backup)

                    if workers_manager.args.debug:
                        for worker in self.workers:
                            if worker.is_online():
                                continue
                                # worker.getDiscoveryServerLog(workers_manager.args.debug)

                    # Report to main shell a script to execute
                    # Delay this to the very end when every worker has received the subgraph
                    print(to_shell(main_graph, workers_manager.args))
                    script_fname = to_shell_file(main_graph, workers_manager.args)
                    log("Master node graph storeid in ", script_fname)
                    response_msg = f"OK {script_fname}"
                    dspash_socket.respond(response_msg, conn)
                except Exception as e:
                    print(e)


            else:
                raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    WorkersManager().run()