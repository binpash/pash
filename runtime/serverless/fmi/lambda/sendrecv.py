import sys

sys.path.append("/opt/python")

import fmi

node_id_str, filename = sys.argv[1:]

node_id = int(node_id_str)

num_nodes = 2
comm = fmi.Communicator(node_id, num_nodes, "fmi.json", "Test", 512)

comm.hint(fmi.hints.fast)

batch = 100000

cnt = 0

if node_id == 0:
    with open(filename, "rb") as f:
        while True:  # Change the loop condition to read until the end of the file
            x = f.read(batch)
            if not x:
                comm.send(batch * [-1], 1, fmi.types(fmi.datatypes.int_list, batch))
                break  # Exit the loop if no more data is read
            send = list(x) + (batch - len(x)) * [-1]
            comm.send(send, 1, fmi.types(fmi.datatypes.int_list, batch))

            cnt += 1
            # print(x)

elif node_id == 1:
    with open(filename, "wb") as f:
        while True:
            x = comm.recv(0, fmi.types(fmi.datatypes.int_list, batch))
            if x[0] == -1:
                break
            write = [item for item in x if item != -1]
            f.write(bytes(write))

            cnt += 1

    print(f"Received {cnt} times")
