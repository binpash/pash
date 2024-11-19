import os
import sys

filename, batch_size = sys.argv[1:3]
batch_size = int(batch_size)

batch = []
# for file in os.listdir(dirname):
#     batch.append(file)
#     if len(batch) == batch_size:
#         print("*".join(batch))
#         batch = []
# if len(batch) > 0:
#     print("*".join(batch))
with open(filename) as f:
    for line in f.readlines():
        batch.append(line.strip())
        if len(batch) == batch_size:
            print("*".join(batch))
            batch = []
if len(batch) > 0:
    print("*".join(batch))