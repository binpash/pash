import socket
import os
import time
import queue
import pickle
import json

from dspash.socket_utils import (
    SocketManager,
    encode_request,
    decode_request,
    send_msg,
    recv_msg,
)
from util import log
from dspash.ir_helper import prepare_graph_for_remote_exec, to_shell_file
from dspash.utils import read_file
import config
import copy

PORT = 65425  # Port to listen on (non-privileged ports are > 1023)


class WorkerConnection:
    def __init__(self, host, port):
        self._host = socket.gethostbyaddr(host)[2][
            0
        ]  # get ip address in case host needs resolving
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

    def send_graph_exec_request(self, graph, shell_vars, functions) -> bool:
        request_dict = {
            "type": "Exec-Graph",
            "graph": graph,
            "functions": functions,
            "shell_variables": None,  # Doesn't seem needed for now
        }
        request = encode_request(request_dict)
        # TODO: do I need to open and close connection?
        send_msg(self._socket, request)
        # TODO wait until the command exec finishes and run this in parallel?
        response_data = recv_msg(self._socket)
        if not response_data or decode_request(response_data)["status"] != "OK":
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


class WorkersManager:
    def __init__(self, workers: WorkerConnection = None):
        self.workers = [] if workers is None else workers
        self.host = socket.gethostbyname(socket.gethostname())
        self.args = copy.copy(config.pash_args)
        # Required to create a correct multi sink graph
        self.args.termination = ""

    def get_worker(self, fids=None) -> WorkerConnection:
        if not fids:
            fids = []

        best_worker = None  # Online worker with least work
        for worker in self.workers:
            if not worker.is_online():
                continue

            # Skip if any provided fid isn't available on the worker machine
            if any(map(lambda fid: not fid.is_available_on(worker.host()), fids)):
                continue

            if (
                best_worker is None
                or best_worker.get_running_processes() > worker.get_running_processes()
            ):
                best_worker = worker

        if best_worker == None:
            raise Exception("no workers online where the date is stored")

        return best_worker

    def add_worker(self, host, port):
        self.workers.append(WorkerConnection(host, port))

    def add_workers_from_cluster_config(self, config_path):
        with open(config_path, "r") as f:
            cluster_config = json.load(f)

        workers = cluster_config["workers"].values()
        for worker in workers:
            host = worker["host"]
            port = worker["port"]
            self.add_worker(host, port)

    def run(self):
        workers_manager = self
        workers_manager.add_workers_from_cluster_config(
            os.path.join(config.PASH_TOP, "cluster.json")
        )

        dspash_socket = SocketManager(os.getenv("DSPASH_SOCKET"))
        while True:
            request, conn = dspash_socket.get_next_cmd()
            if request.startswith("Done"):
                dspash_socket.close()
                break
            elif request.startswith("Exec-Graph"):
                args = request.split(":", 1)[1].strip()
                filename, declared_functions_file = args.split()

                (
                    worker_subgraph_pairs,
                    shell_vars,
                    main_graph,
                ) = prepare_graph_for_remote_exec(filename, self.get_worker)
                script_fname = to_shell_file(main_graph, self.args)
                log("Master node graph stored in ", script_fname)

                # Read functions
                log("Functions stored in ", declared_functions_file)
                declared_functions = read_file(declared_functions_file)

                # Report to main shell a script to execute
                response_msg = f"OK {script_fname}"
                dspash_socket.respond(response_msg, conn)

                # Execute subgraphs on workers
                for worker, subgraph in worker_subgraph_pairs:
                    worker.send_graph_exec_request(
                        subgraph, shell_vars, declared_functions
                    )
            else:
                raise Exception(f"Unknown request: {request}")


if __name__ == "__main__":
    WorkersManager().run()
