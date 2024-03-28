import socket
import os
import time
import queue
import pickle
import json
from uuid import uuid4

from dspash.socket_utils import SocketManager, encode_request, decode_request, send_msg, recv_msg
from util import log
from dspash.ir_helper import prepare_graph_for_remote_exec, to_shell_file, add_debug_flags, is_merger_subgraph, has_merger_subgraph, get_main_writer_graphs
from dspash.utils import read_file
from dspash.hdfs_utils import start_hdfs_deamon, stop_hdfs_deamon
import config 
import copy
import requests
from ir_to_ast import to_shell
import time
import collections

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

    def send_graph_exec_request(self, graph, shell_vars, functions, args) -> bool:
        request_dict = { 'type': 'Exec-Graph',
                        'graph': graph,
                        'functions': functions,
                        'shell_variables': None, # Doesn't seem needed for now
                        'debug': None,
                        'kill_target': None,
                    }
        if args.debug:
            request_dict['debug'] = {'name': self.name, 'url': f'{DEBUG_URL}/putlog'}

        if args.kill:
            request_dict['kill_target'] = socket.gethostbyaddr(args.kill)[2][0]
            log(f"Kill target: {request_dict['kill_target']}")

        request = encode_request(request_dict)
        #TODO: do I need to open and close connection?
        send_msg(self._socket, request)
        # TODO wait until the command exec finishes and run this in parallel?
        response_data = recv_msg(self._socket)
        if not response_data or decode_request(response_data)['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
        else:
            # self._running_processes += 1 #TODO: decrease in case of failure or process ended
            response = decode_request(response_data)
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
    def __init__(self, workers: WorkerConnection = []):
        self.workers = workers
        self.host = socket.gethostbyname(socket.gethostname())
        self.args = copy.copy(config.pash_args)
        # Required to create a correct multi sink graph
        self.args.termination = ""
        self.request_args = [] 
        self.dspash_socket = None
        # metadata
        self.subgraph_id_to_subgraph = {}
        self.worker_to_subgraphs_id = collections.defaultdict(list)
        self.subgraph_id_to_worker = collections.defaultdict(WorkerConnection)

    def init_metadata(self, worker_subgraph_pairs):
        # Excluding main_graph for now
        for worker, subgraph in worker_subgraph_pairs:
            subgraph_id = subgraph.id
            self.subgraph_id_to_subgraph[subgraph_id] = subgraph
            self.worker_to_subgraphs_id[worker].append(subgraph_id)
            self.subgraph_id_to_worker[subgraph_id] = worker

    def get_subgraphs_for_worker(self, worker: WorkerConnection):
        return [self.subgraph_id_to_subgraph[subgraph_id] for subgraph_id in self.worker_to_subgraphs_id[worker]]


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

        addrs = {conn.host() for conn in self.workers}
        start_hdfs_deamon(10, addrs, self.addr_added, self.addr_removed)

    def addr_added(self, addr: str):
        log(f"Added {addr} to active nodes")

    # This function is called when a datanode is found to be dead. 
    # Fault tolerance mechanism is triggered
    def addr_removed(self, addr: str):
        log(f"Removed {addr} from active nodes")
        worker = self.get_worker_by_host(addr.split(':')[0])
        worker.close()
        log(f"{worker} closed")
        # Check if the crashed worker contains a merger subgraph
        merger_node_crash = has_merger_subgraph(self.get_subgraphs_for_worker(worker))
        if merger_node_crash:
            self.schedule_graphs(True, False)
        else:
            self.schedule_graphs(True, True)
        log(f"Fault tolerance mechanism triggered")

    def get_worker_by_host(self, host):
        for worker in self.workers:
            if worker.host() == host:
                return worker
        return None
    
    # If skip_merger_subgraph, it means 1) we don't schedule the merger subgraph, 2) no need to re-execute part of main_graph
    def schedule_graphs(self, executingFT: bool, skip_merger_subgraph, conn: socket.socket=None):
        filename, declared_functions_file = self.request_args.split()
        reExecuting = False
        workers_manager = self
        try:
            # In the naive fault tolerance, if a crash happens, 
            # we'll re-split the IR and do it again until scheduling is done without any crash.
            worker_subgraph_pairs, shell_vars, main_graph = prepare_graph_for_remote_exec(filename, workers_manager.get_worker, self.request_hash)
            log(f"num subgraphs: {len(worker_subgraph_pairs)}")

            self.init_metadata(worker_subgraph_pairs)
            # Read functions
            log("Functions stored in ", declared_functions_file)
            declared_functions = read_file(declared_functions_file)

            if not executingFT:
                add_debug_flags(main_graph)

                # Only send main_graph to client once. Re-executing logic will be handled in datastream.go
                script_fname = to_shell_file(main_graph, workers_manager.args)
                log("Master node graph stored in ", script_fname)
            
                # Report to main shell a script to execute
                log(to_shell(main_graph, workers_manager.args))
                log("-------------------------------------")
                response_msg = f"OK {script_fname}"
                self.dspash_socket.respond(response_msg, conn)
                log("-------------------------------------")
            else:
                if not skip_merger_subgraph:
                    # if merger main_graph is crashed, we need to re-execute part of main_graph on client's node

                    # when a main_graph is a "merger subgraph", it doesn' mean it has a merger command
                    #       it means in addition to reading from merger subgraph with remote_read, it 
                    #       also sends some local data (e.g. from local file on client's node) to the [actual] merger subgraph

                    # Check if there's a part of main_graph that also needs to be re-executed
                    # e.g. if main_graph contains 1) datastream write < cat dict.txt, 2) datastream read, the first part 1) needs to be re-executed
                    
                    # think of this as "customer service" after sending main_graph to the client initially
                    # start a new connection to the socket, send a graph for re-execution when fault happens
                    #                 and there's a need to re-execute part of main_grpah on the client node
                    main_writer_graphs = get_main_writer_graphs(main_graph)
                    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
                        sock.connect(os.getenv('DSPASH_REEXEC_SOCKET'))
                        script_fnames = []
                        for main_writer_graph in main_writer_graphs:
                            if not main_writer_graph.empty():
                                log("main_writer_graph")
                                log(to_shell(main_writer_graph, workers_manager.args))
                                script_fname = to_shell_file(main_writer_graph, workers_manager.args)
                                script_fnames.append(script_fname)
                        try:
                            script_fnames_str = " ".join(script_fnames)
                            message = f"REEXEC {script_fnames_str}\n"
                            bytes_message = message.encode('utf-8')
                            sock.sendall(bytes_message)
                        except Exception as e:
                            log(e)
                        log("-------------------------------------")
            # Execute subgraphs on workers
            for worker, subgraph in worker_subgraph_pairs: 
                if skip_merger_subgraph and is_merger_subgraph(subgraph):
                    log("Skipping merger subgraph!")
                    continue                           
                try:
                    log(worker.name, worker)
                    log(to_shell(subgraph, workers_manager.args))
                    log("-------------------------------------")
                    worker.send_graph_exec_request(subgraph, shell_vars, declared_functions, workers_manager.args)
                except Exception as e:
                    # worker timeout
                    worker.close()
                    log(f"{worker} closed")
                    break
        except Exception as e:
            print(e)


    def run(self):
        workers_manager = self
        workers_manager.add_workers_from_cluster_config(os.path.join(config.PASH_TOP, 'cluster.json'))
        if workers_manager.args.debug:
            try:
                requests.post(f'{DEBUG_URL}/clearall') # clears all the debug server logs
            except Exception as e:
                log(f"Failed to connect to debug server with error {e}\n")
                workers_manager.args.debug = False # Turn off debugging

        self.dspash_socket = SocketManager(os.getenv('DSPASH_SOCKET'))
        while True:
            request, conn = self.dspash_socket.get_next_cmd()
            if request.startswith("Done"):
                self.dspash_socket.close()
                stop_hdfs_deamon()
                break
            elif request.startswith("Exec-Graph"):
                self.request_args = request.split(':', 1)[1].strip()
                self.request_hash = str(uuid4())
                self.schedule_graphs(False, False, conn)
            else:
                stop_hdfs_deamon()
                raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    WorkersManager().run()