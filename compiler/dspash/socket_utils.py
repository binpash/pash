import socket
import os
import json
import shlex
import subprocess
import pickle
import struct


def send_msg(sock, msg):
    # Prefix each message with a 4-byte length (network byte order)
    msg = struct.pack(">I", len(msg)) + msg
    sock.sendall(msg)


def recv_msg(sock):
    # Read message length and unpack it into an integer
    raw_msglen = recvall(sock, 4)
    if not raw_msglen:
        return None
    msglen = struct.unpack(">I", raw_msglen)[0]
    # Read the message data
    return recvall(sock, msglen)


def recvall(sock, n):
    # Helper function to recv n bytes or return None if EOF is hit
    data = bytearray()
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return None
        data.extend(packet)
    return data


def encode_request(obj: dict):
    return pickle.dumps(obj)


def decode_request(b: bytes):
    return pickle.loads(b)


## TODO: SocketManager might need to handle errors more gracefully
class SocketManager:
    def __init__(self, server_address):
        ## Configure them outside
        self.buf_size = 8192

        # Make sure the socket does not already exist
        ## TODO: Is this necessary?
        try:
            os.unlink(server_address)
        except OSError:
            if os.path.exists(server_address):
                raise
        # log("SocketManager: Made sure that socket does not exist")

        # Create a UDS socket
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        # log("SocketManager: Created socket")

        self.sock.bind(server_address)
        # log("SocketManager: Successfully bound to socket")

        ## TODO: Check if we need to configure the back# log
        self.sock.listen()
        # log("SocketManager: Listenting on socket")

    def get_next_cmd(self):
        connection, client_address = self.sock.accept()
        data = connection.recv(self.buf_size)

        ## TODO: This could be avoided for efficiency
        str_data = data.decode("utf-8")
        # log("Received data:", str_data)
        ## TODO: Lift this requirement if needed
        ##
        ## We need to ensure that we read a command at once or the command was empty (only relevant in the first invocation)
        assert str_data.endswith("\n") or str_data == ""

        return str_data, connection

    ## This method respond to the connection we last got input from
    ## In the case of the UnixPipes, we don't have any state management here
    ##   since all reads/writes go to/from the same fifos
    def respond(self, message, connection):
        bytes_message = message.encode("utf-8")
        connection.sendall(bytes_message)
        connection.close()

    def close(self):
        self.sock.close()
        # log("SocketManager: Closed")
