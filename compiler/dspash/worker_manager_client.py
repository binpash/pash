from socket_utils import SocketManager

import socket
import sys
import subprocess
import threading

def log(*args, end='\n', level=1):
    ## If the debug logging level is at least
    ## as high as this log message.
    print("worker_manager_client: ", *args, file=sys.stderr, end=end, flush=True)


# Send requestto worker_manager
def send_request(ir_file, declared_functions, dspash_socket, dspash_reexec_socket):
    try:
        # Listen for responses in the main thread
        # listen_thread = threading.Thread(target=listen_responses, args=(dspash_reexec_socket, ))
        listen_thread = threading.Thread(target=listen_responses, args=(dspash_reexec_socket, ), daemon=True)
        listen_thread.start()
        # Connect to the worker manager via Unix socket
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
            sock.connect(dspash_socket)
            # Send request to worker manager
            message = f"Exec-Graph: {ir_file} {declared_functions}\n"
            sock.sendall(message.encode('utf-8'))
            # Receive response from worker manager
            response = sock.recv(4096).decode('utf-8').split(' ')
            status, script_to_execute = response[0], response[1]
            cmd = cmd = f"source {script_to_execute}"
            rc = subprocess.Popen(
                cmd, executable="/bin/bash", shell=True)
            rc.wait()
            
    except Exception as e:
        log(f"Error: {e}")

# Receive response (main_graph) from worker_manager
# Including potential re-execution of main_graph
def listen_responses(dspash_reexec_socket):
    sock = SocketManager(dspash_reexec_socket)
    reexec_rcs = []
    while True:
        request, conn = sock.get_next_cmd()
        req_type, scripts_to_execute = request.split(' ')[0], request.split(' ')[1]
        if req_type == "REEXEC":
            for script_to_execute in scripts_to_execute.split(' '):
                cmd = f"source {script_to_execute}"
                rc = subprocess.Popen(cmd, executable="/bin/bash", shell=True)
                reexec_rcs.append(rc)
            break
        else:
            log("Error:", script_to_execute)
    log("Received all re-execution requests, wait for them to finish.")
    for rc in reexec_rcs:
        rc.wait()
    log("All re-execution commands finished!")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python3 worker_manager_client.py ir_file declared_functions socket_path reexec-socket_path")
        sys.exit(1)

    # Assign command-line arguments to variables
    ir_file = sys.argv[1]
    declared_functions = sys.argv[2]
    dspash_socket = sys.argv[3]
    dspash_reexec_socket = sys.argv[4]

    send_request(ir_file, declared_functions, dspash_socket, dspash_reexec_socket)
