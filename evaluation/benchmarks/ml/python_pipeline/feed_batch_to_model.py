import sys
from torch import load, save

model = load(sys.argv[1])
batch = load(sys.argv[2])
batch = batch.view(batch.shape[0], -1)

output = model(batch)
save(output, 'output.pt')
print('output.pt')