import socket
import os
import time
import queue
import pickle

from dspash.socket_utils import SocketManager, encode_request, decode_request, send_msg, recv_msg
from util import log
from dspash.ir_helper import prepare_graph_for_remote_exec
import config 

PORT = 65425        # Port to listen on (non-privileged ports are > 1023)
START_PORT = 57000

class WorkerConnection:
    def __init__(self, host, port):
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

    def send_graph_exec_request(self, graph, shell_vars) -> bool:
        request_dict = { 'type': 'Exec-Graph',
                        'graph': graph,   
                        'shell_variables': None # Doesn't seem needed for now     
                    }
        request = encode_request(request_dict)
        #TODO: do I need to open and close connection?
        send_msg(self._socket, request)
        # TODO wait until the command exec finishes and run this in parallel?
        response_data = recv_msg(self._socket)
        if not response_data or decode_request(response_data)['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
        else:
            self._running_processes += 1 #TODO: decrease in case of failure or process ended
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
        
        config.next_available_port = START_PORT # TODO: improve to per worker config

    def get_worker(self, fids) -> WorkerConnection:
        best_worker = None  # Online worker with least work
        for worker in self.workers:
            if not worker.is_online():
                continue
            # Can't use this worker if data is not availabe on it
            for fid in fids:
                # TODO: Can this be optimized by only checking sink and source fids?
                if not fid.is_available_on(worker.host()):
                    continue

            if best_worker is None or best_worker.get_running_processes() > worker.get_running_processes():
                best_worker = worker

        if best_worker == None:
            raise Exception("no workers online where the date is stored")

        return best_worker

    def add_worker(self, host, port):
        self.workers.append(WorkerConnection(host, port))

    def run(self):
        workers_manager = self
        workers_manager.add_worker(self.host, 65432)
        workers_manager.add_worker(self.host, 65431)

        dspash_socket = SocketManager(os.getenv('DSPASH_SOCKET'))
        while True:
            request, conn = dspash_socket.get_next_cmd()
            if request.startswith("Done"):
                dspash_socket.close()
                break
            elif request.startswith("Exec-Graph"):
                filename = request.split(':', 1)[1].strip()
                graphs, shell_vars, final_graph_file = prepare_graph_for_remote_exec(filename)
                log("Master node graph stored in ", final_graph_file)
                # Report to main shell where to read from
                response_msg = f"OK {final_graph_file}"
                dspash_socket.respond(response_msg, conn) 
                # Execute subgraphs on workers
                for i in range(len(graphs) - 1, -1, -1):
                    subgraph = graphs[i]
                    fids = subgraph.all_fids()
                    worker = workers_manager.get_worker(fids)
                    worker.send_graph_exec_request(subgraph, shell_vars)
            else:
                raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    WorkersManager().run()
