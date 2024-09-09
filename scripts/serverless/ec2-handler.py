import socket
import threading
import subprocess
import uuid
import os
import json

# Function to handle each client connection
def handle_client(client_socket, client_address, tmp_dir):
    print(f"[EC2 Handler] Accepted connection from {client_address}")
    try:
        request = b""
        while True:
            chunk = client_socket.recv(1024)
            if not chunk:
                break
            request += chunk
        if not request:
            # print(f"[EC2 Handler] Closing connection with {client_address}")
            client_socket.close()
            return
        event = json.loads(request.decode("utf-8"))
        data = event["data"]
        id_ = event["id"]
        scripts_dict = json.loads(data)
        with open(f"{tmp_dir}/data-{id_}", "w") as f:
            f.write(data)
        print("[EC2 Handler] Executing script ID", id_)
        # print(scripts_dict[id_])
        with open(f"{tmp_dir}/script-{id_}.sh", "w") as f:
            f.write(scripts_dict[id_])
        process = subprocess.run(
            ["/bin/bash", f"{tmp_dir}/script-{id_}.sh", f"{tmp_dir}/data-{id_}"]
        )
        print(f"[EC2 Handler] script {id_} execution return code: {process.returncode}")
    except Exception as e:
        print(f"[EC2 Handler] Error handling client {client_address}: {e}")
    finally:
        print(f"[EC2 Handler] Closing connection with {client_address}")
        client_socket.close()

# Main function to start the server
def start_server(host='0.0.0.0', port=9999):
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((host, port))
    server.listen(5)  # Listen for incoming connections (5 is the max number of queued connections)
    uid = str(uuid.uuid4())
    tmp_dir = os.path.dirname(f"/tmp/pash_{uid}/")
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    print(f"[EC2 Handler] Server listening on {host}:{port}, tmp directory: {tmp_dir}")

    while True:
        client_socket, client_address = server.accept()  # Accept a new connection
        client_handler = threading.Thread(target=handle_client, args=(client_socket, client_address, tmp_dir))
        client_handler.start()  # Start a new thread to handle the client connection

if __name__ == "__main__":
    start_server(port=9999)