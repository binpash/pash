from collections import defaultdict
import socket
import os
from threading import Event, Thread
import time
import queue
import pickle
import json
from typing import List, Tuple
from uuid import UUID

from dspash.socket_utils import SocketManager, encode_request, decode_request, send_msg, recv_msg
from ir import IR
from util import log
from ir_to_ast import to_shell
from dspash.ir_helper import prepare_graph_for_remote_exec, to_shell_file, add_debug_flags, get_main_writer_graphs
from dspash.utils import read_file
from dspash.hdfs_utils import start_hdfs_daemon, stop_hdfs_daemon
import config 
import copy
import requests

HOST = socket.gethostbyname(socket.gethostname())
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

    def send_kill_all_subgraphs_request(self) -> bool:
        request_dict = { 'type': 'Kill-All-Subgraphs' }
        request = encode_request(request_dict)
        send_msg(self._socket, request)
        response_data = recv_msg(self._socket)
        if not response_data or decode_request(response_data)['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
        else:
            # response = decode_request(response_data)
            return True

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

        log(f"Sending graph {graph.id} to {self.name}")
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
        self._socket.send("Done")
        self._socket.close()

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

        self.worker_subgraph_pairs: List[Tuple[WorkerConnection, IR]] = []
        self.main_graph = IR({}, {})
        self.shell_vars = None
        self.uuid_to_graphs = {}
        self.graph_to_uuid = defaultdict(list)
        self.declared_functions = None

        self.daemon_quit = Event()
        self.uuid_to_read_graph_id = {}
        self.uuid_to_write_graph_id = {}
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.bind((HOST, PORT))
        self.s.listen()
        log(f"Worker manager on {HOST}:{PORT}")
        Thread(target=self.__daemon, daemon=True).start()

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
        start_hdfs_daemon(10, addrs, self.addr_added, self.addr_removed)

    def addr_added(self, addr: str):
        log(f"Added {addr} to active nodes")
        for worker in self.workers:
            if worker.host() == addr:
                worker._online = True

    def addr_removed(self, addr: str):
        log(f"Removed {addr} from active nodes")
        for worker in self.workers:
            if worker.host() == addr:
                worker._online = False

        # merger node crash
        for worker, subgraph in self.worker_subgraph_pairs:
            if subgraph.merger and worker.host() == addr:
                self.handle_merger_crash(addr)
                return

        # regular node crash
        self.handle_regular_crash(addr)

    def handle_regular_crash(self, addr: str):
        log("Regular node crashed, handling it, ft mode is", self.args.ft)
        try:
            if self.args.ft == "base" or self.args.ft == "optimized":
                # iterate over the copy of the list to avoid modifying it while iterating
                for worker, subgraph in self.worker_subgraph_pairs[:]:
                    # second condition is to avoid sending the subgraph 
                    # if it's already sent all of it's outputs to another worker
                    if worker.host() == addr and self.graph_to_uuid[subgraph.id]:
                        self.worker_subgraph_pairs.remove((worker, subgraph))
                        subgraph_critical_fids = list(filter(lambda fid: fid.has_remote_file_resource(), subgraph.all_fids()))
                        new_worker = self.get_worker(subgraph_critical_fids)
                        new_worker._running_processes += 1
                        self.worker_subgraph_pairs.append((new_worker, subgraph))

                        new_worker.send_graph_exec_request(subgraph, self.shell_vars, self.declared_functions, self.args)
                        log(f"Sent subgraph {subgraph.id} to {new_worker}")
            
            if self.args.ft == "naive":
                # Regular crash isn't different from merger crash with naive ft mode
                self.handle_merger_crash(addr)
            
        except Exception as e:
            log(f"Failed to handle regular crash with error {e}")

    def handle_merger_crash(self, addr):
        log("Merger node crashed, handling it, ft mode is", self.args.ft)
        try:
            if self.args.ft == "base" or self.args.ft == "naive":
                log(f"FT mode is {self.args.ft}, killing all subgraphs")
                for worker in self.workers:
                    worker: WorkerConnection
                    worker.send_kill_all_subgraphs_request()
                    log(f"Sent kill all subgraphs request to {worker}")
                log("Killed all subgraphs, will send all subgraphs again")

                self.generate_or_reset_mapping()

                # iterate over the copy of the list to avoid modifying it while iterating
                for worker, subgraph in self.worker_subgraph_pairs[:]:
                    if worker.host() == addr:
                        self.worker_subgraph_pairs.remove((worker, subgraph))
                        subgraph_critical_fids = list(filter(lambda fid: fid.has_remote_file_resource(), subgraph.all_fids()))
                        new_worker = self.get_worker(subgraph_critical_fids)
                        new_worker._running_processes += 1
                        self.worker_subgraph_pairs.append((new_worker, subgraph))
                        worker = new_worker

                    worker.send_graph_exec_request(subgraph, self.shell_vars, self.declared_functions, self.args)
                log("Sent all subgraphs again")

            if self.args.ft == "optimized":
                pass

            # Check if there's a part of main_graph that also needs to be re-executed
            # e.g. if main_graph contains 1) datastream write < cat dict.txt, 2) datastream read, the first part 1) needs to be re-executed

            # think of this as "customer service" after sending main_graph to the client initially
            # start a new connection to the socket, send a graph for re-execution when fault happens
            #                 and there's a need to re-execute part of main_grpah on the client node
            main_writer_graphs = get_main_writer_graphs(self.main_graph)
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
                sock.connect(os.getenv('DSPASH_REEXEC_SOCKET'))
                script_fnames = []
                for main_writer_graph in main_writer_graphs:
                    if not main_writer_graph.empty():
                        log("main_writer_graph")
                        log(to_shell(main_writer_graph, self.args))
                        script_fname = to_shell_file(main_writer_graph, self.args)
                        script_fnames.append(script_fname)
                if script_fnames:
                    try:
                        script_fnames_str = " ".join(script_fnames)
                        message = f"REEXEC {script_fnames_str}\n"
                        bytes_message = message.encode('utf-8')
                        sock.sendall(bytes_message)
                    except Exception as e:
                        log(e)
                    log("-------------------------------------")

        except Exception as e:
            log(f"Failed to handle merger crash with error {e}")

    def __daemon(self):
        self.s.settimeout(1)  # Set a timeout of 1 second
        while not self.daemon_quit.is_set():
            try:
                conn, addr = self.s.accept()
                Thread(target=self.__manage_connection, args=[conn, addr]).start()
            except socket.timeout:
                continue

    def __manage_connection(self, conn: socket, addr):
        # 1 byte for read or write and 16 byte for uuid
        data = conn.recv(17)
        # ideally do this in a loop, but it's extemely unlikely that the data will be split
        assert len(data) == 17
        # Read the first byte to get if the request is from read or write client
        read_client = True if data[0] == 0 else False

        uuid = UUID(bytes=data[1:])
        if read_client:
            responsible_graph = self.uuid_to_graphs[uuid][0]
            self.graph_to_uuid[responsible_graph].remove(uuid)
            log(f"Removed {uuid} from graph_to_uuid[{responsible_graph}]")

    def generate_or_reset_mapping(self):
        self.graph_to_uuid.clear()
        # These graphs are responsible for producing these outputs
        for uuid, (from_graph, _) in self.uuid_to_graphs.items():
            self.graph_to_uuid[from_graph].append(uuid)

    def run(self):
        self.add_workers_from_cluster_config(os.path.join(config.PASH_TOP, 'cluster.json'))
        if self.args.debug:
            try:
                requests.post(f'{DEBUG_URL}/clearall') # clears all the debug server logs
            except Exception as e:
                log(f"Failed to connect to debug server with error {e}\n")
                self.args.debug = False # Turn off debugging

        dspash_socket = SocketManager(os.getenv('DSPASH_SOCKET'))
        while True:
            request, conn = dspash_socket.get_next_cmd()
            if request.startswith("Done"):
                dspash_socket.close()
                stop_hdfs_daemon()
                self.daemon_quit.set()
                break
            elif request.startswith("Exec-Graph"):
                args = request.split(':', 1)[1].strip()
                filename, declared_functions_file = args.split()

                self.worker_subgraph_pairs, self.shell_vars, self.main_graph, self.uuid_to_graphs = prepare_graph_for_remote_exec(filename, self.get_worker)
                script_fname = to_shell_file(self.main_graph, self.args)
                if self.args.debug:
                    add_debug_flags(self.main_graph)
                log("Master node graph stored in ", script_fname)

                self.generate_or_reset_mapping()

                # Read functions
                log("Functions stored in ", declared_functions_file)
                self.declared_functions = read_file(declared_functions_file)

                # Report to main shell a script to execute
                log("main graph")
                log(to_shell(self.main_graph, self.args))
                log("-------------------------------------")
                response_msg = f"OK {script_fname}"
                dspash_socket.respond(response_msg, conn)

                # Execute subgraphs on workers
                for worker, subgraph in self.worker_subgraph_pairs:
                    worker: WorkerConnection
                    worker.send_graph_exec_request(subgraph, self.shell_vars, self.declared_functions, self.args)
            else:
                stop_hdfs_daemon()
                self.daemon_quit.set()
                raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    WorkersManager().run()
