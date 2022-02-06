import socket
from threading import Thread
from socket_utils import SocketManager, encode_request, decode_request
import os
import time
HOST = '0.0.0.0'
PORT = 65425        # Port to listen on (non-privileged ports are > 1023)

class WorkerConnection:
    def __init__(self, host, port):
        self.host = host
        self.port = port
        self.running_processes = 0
        self.online = True
        # assume client service is running, can add a way to activate later
        try:
            self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.socket.connect((host, port))
        except:
            self.online = False

    def is_online(self):
        # TODO: create a ping to confirm is online
        return self.online

    def get_running_processes(self):
        # request_dict = { 'type': 'query',
        #                 'fields': ['running_processes']
        # }
        # with self.socket.connect((self.host, self.port)):
        #     self.socket.send(request)
        #     # TODO wait until the command exec finishes and run this in parallel?
        #     answer = self.socket.recv(1024)
        return self.running_processes
            
    def send_exec_request(self, cmd) -> (bool, int, int):
        request_dict = { 'type': 'exec',
                    'cmd': cmd,
                   }
        request = encode_request(request_dict)
        #TODO: do I need to open and close connection?
        self.socket.send(request)
        # TODO wait until the command exec finishes and run this in parallel?
        response_data = self.socket.recv(4096)
        if not response_data or decode_request(response_data)['status'] != "OK":
            raise Exception(f"didn't recieved ack on request {response_data}")
            return False, -1, -1
        else:
            self.running_processes += 1 #TODO: decrease in case of failure or process ended
            response = decode_request(response_data)
            return True, response['body']['from_port'], response['body']['to_port']
            

    def close(self):
        self.socket.send("Done")
        self.socket.close()

    def _wait_ack(self):
        confirmation = self.socket.recv(4096)
        if not confirmation or decode_request(confirmation).status != "OK":
            return False
            # log(f"Confirmation not recieved {confirmation}")
        else:
            return True

    def __str__(self):
        return f"Worker {self.host}:{self.port}"


class WorkersManager():
    def __init__(self, workers: WorkerConnection = []):
        self.workers = workers

    def get_worker(self) -> WorkerConnection:
        best_worker = None  # Online worker with least work
        for worker in self.workers:
            if not worker.is_online():
                continue
            if best_worker is None or best_worker.running_processes > worker.running_processes:
                best_worker = worker

        if best_worker == None:
            raise Exception("no workers online")

        return best_worker

    def add_worker(self, host, port):
        self.workers.append(WorkerConnection(host, port))

def main():
    workers_manager = WorkersManager()
    workers_manager.add_worker('0.0.0.0', 65432)
    workers_manager.add_worker('0.0.0.0', 65431)
    request1 = "Get-Worker: $PASH_TOP/runtime/r_wrap grep 'Bell'"
    request2 = "Get-Worker: $PASH_TOP/runtime/r_wrap grep 'Bell'"
    dspash_socket = SocketManager(os.getenv('DSPASH_SOCKET'))
    while True:
        request, conn = dspash_socket.get_next_cmd()
        if request.startswith("Done"):
            dspash_socket.close()
            break
        elif request.startswith("Get-Worker"):
            cmd = request.split(':', 1)[1].strip()
            worker = workers_manager.get_worker()
            ok, from_port, to_port = worker.send_exec_request(cmd)
            # time.sleep(5)
            if not ok: 
                raise Exception("Failed sending exec request")
            #TODO: send port
            response_msg = f"OK {worker.host} {from_port} {to_port}"
            dspash_socket.respond(response_msg, conn)
        else:
            raise Exception(f"Unknown request: {request}")
        
if __name__ == "__main__":
    main()
