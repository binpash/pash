import os
import socket

import config
from util import log


def success_response(string):
    return f"OK: {string}\n"


def error_response(string):
    return f"ERROR: {string}\n"


class UnixPipeReader:
    def __init__(self, in_filename, out_filename, blocking=True):
        self.in_filename = in_filename
        self.out_filename = out_filename
        self.buffer = ""
        self.blocking = blocking
        if not self.blocking:
            # Non blocking mode shouldn't be used in production. It's only used experimentally.
            log("Reader initialized in non-blocking mode")
            self.fin = open(self.in_filename)
        else:
            log("Reader initialized in blocking mode")

    ## This is necessary here to ensure that get_next_cmd is blocking een though the underlying API is non-blocking.
    def get_next_cmd(self):
        cmd = ""
        ## TODO: Remove the non-blocking control flow since it doesn't make sense to busy wait
        if not self.blocking:
            while not cmd:
                cmd = self.get_next_cmd_aux()
        else:
            cmd = self.get_next_cmd_aux()
        return cmd

    def get_next_cmd_aux(self):
        """
        This method return depends on the reading mode. In blocking mode this method will
        return the next full command and if there is no command it will wait until a full command is recieved.
        In non blocking mode it would either a full command or an empty string if a full command isn't available yet.
        This command keeps a state of the remaining data which is used in each subsequent call to this method.
        """
        input_buffer = ""
        if self.buffer:
            # Don't wait on fin if cmd buffer isn't empty
            log(
                "Reader buffer isn't empty. Using it instead of reading new data for the next command"
            )
            input_buffer = self.buffer
        else:
            log("Reader buffer is empty. Reading new data from input fifo")
            if self.blocking:
                with open(self.in_filename) as fin:
                    # This seems to be necessary for reading the full data.
                    # It seems like slower/smaller machines might not read the full data in one read
                    while True:
                        data = fin.read()
                        if len(data) == 0:
                            break
                        input_buffer += data
            else:
                input_buffer = self.fin.read()

        log("Input buffer:", input_buffer)
        if "\n" in input_buffer:
            cmd, rest = input_buffer.split("\n", 1)  # split on the first \n only
            self.buffer = rest
        else:
            cmd = input_buffer
            self.buffer = ""

        cmd = cmd.rstrip()
        log("Reader returned cmd:", cmd)
        return cmd

    ## This method respond to the connection we last got input from
    ## In the case of the UnixPipes, we don't have any state management here
    ##   since all reads/writes go to/from the same fifos
    def respond(self, message):
        fout = open(self.out_filename, "w")
        fout.write(message)
        fout.flush()
        fout.close()

    ## This method doesn't do anything for unix pipe reader since we always read and write
    ## to and from the same fifos
    def close_last_connection(self):
        pass

    def close(self):
        log("Reader closed")
        if not self.blocking:
            self.fin.close()


def unix_socket_send_and_forget(socket_file: str, msg: str):
    try:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(socket_file)
        msg_with_newline = msg + "\n"
        byte_msg = msg_with_newline.encode("utf-8")
        sock.sendall(byte_msg)
        data = sock.recv(config.SOCKET_BUF_SIZE)
        str_data = data.decode("utf-8")
        ## There should be no response on these messages
        assert len(str_data) == 0
    finally:
        log("Sent message:", msg, "to server.", level=1)
        sock.close()


## TODO: Instead of this, think of using a standard SocketServer
##   see: https://docs.python.org/3/library/socketserver.html#module-socketserver
##
## TODO: SocketManager might need to handle errors more gracefully
class SocketManager:
    def __init__(self, socket_addr: str):
        ## Configure them outside
        server_address = socket_addr
        self.buf_size = config.SOCKET_BUF_SIZE

        # Make sure the socket does not already exist
        ## TODO: Is this necessary?
        try:
            os.unlink(server_address)
        except OSError:
            if os.path.exists(server_address):
                raise
        log("SocketManager: Made sure that socket does not exist")

        # Create a UDS socket
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        log("SocketManager: Created socket")

        self.sock.bind(server_address)
        log("SocketManager: Successfully bound to socket")

        ## TODO: Check if we need to configure the backlog
        self.sock.listen()
        log("SocketManager: Listenting on socket")

        ## Connection stack
        self.connections = []

    def get_next_cmd(self):
        connection, client_address = self.sock.accept()
        data = connection.recv(self.buf_size)

        ## TODO: This could be avoided for efficiency
        str_data = data.decode("utf-8")
        log("Received data:", str_data)
        ## TODO: Lift this requirement if needed
        ##
        ## We need to ensure that we read a command at once or the command was empty (only relevant in the first invocation)
        assert str_data.endswith("\n") or str_data == ""

        self.connections.append(connection)
        return str_data

    ## This method respond to the connection we last got input from
    ## In the case of the UnixPipes, we don't have any state management here
    ##   since all reads/writes go to/from the same fifos
    def respond(self, message):
        bytes_message = message.encode("utf-8")
        self.connections[-1].sendall(bytes_message)
        self.close_last_connection()

    ## This method doesn't do anything for unix pipe reader since we always read and write
    ## to and from the same fifos
    def close_last_connection(self):
        # Clean up the connection
        last_connection = self.connections.pop()
        last_connection.close()

    def close(self):
        self.sock.close()
        log("SocketManager: Closed")
