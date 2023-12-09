import sys

sys.path.append("/opt/python")

import fmi

key, fifo, node_id = sys.argv[1:]

node_id = int(node_id)

num_nodes = 2
comm = fmi.Communicator(node_id, num_nodes, "aws/fmi.json", key, 512)

comm.hint(fmi.hints.fast)

batch = 1000000

cnt = 0

def try_send(data, comm):
    while True:
        try:
            comm.send(data, 1, fmi.types(fmi.datatypes.int_list, batch))
            break
        except Exception as e:
            print(f"Failed to send, trying re-sending: {e}", flush=True)

def try_receive(comm):
    while True:
        try:
            data = comm.recv(0, fmi.types(fmi.datatypes.int_list, batch))
            return data
        except Exception as e:
            print(f"Failed to receive, trying re-receiving: {e}", flush=True)

if node_id == 0:
    # dummy send
    try_send([1]*batch, comm)
    with open(fifo, "rb") as f:
        while True:  # Change the loop condition to read until the end of the file
            x = f.read(batch)
            if not x:
                data = batch * [-1]
                print(f"{key}: try sending data with length {len(x)} (all negative), counter {cnt}", flush=True)
                try_send(data, comm)
                print(f"{key}: sent data with length {len(x)} (all negative), counter {cnt}", flush=True)
                break  # Exit the loop if no more data is read
            data = list(x) + (batch - len(x)) * [-1]
            print(f"{key}: try sending data with length {len(x)}, counter {cnt}", flush=True)
            try_send(data, comm)
            print(f"{key}: sent data with length {len(x)}, counter {cnt}", flush=True)
            cnt += 1

elif node_id == 1:
    # dummy recv
    try_receive(comm)
    with open(fifo, "wb") as f:
        while True:
            print(f"{key}: try recving data with counter {cnt}", flush=True)
            x = try_receive(comm)
            if x[0] == -1:
                print(f"{key}: recv data with length {len(x)} (first byte negative), counter {cnt}", flush=True)
                break
            write = [item for item in x if item != -1]
            f.write(bytes(write))
            f.flush()
            print(f"{key}: (flushed) recv data with length {len(write)}, counter {cnt}", flush=True)
            cnt += 1
            # if len(write) < batch:
            #     break
