import sys
from torch import load, save

train_loader = load(sys.argv[1])

i = 0
for batch, labels in train_loader:
    save(batch, f'./intermed/batches/batch_{i}.pt')
    save(labels, f'./intermed/labels/labels_{i}.pt')
    i += 1
print(f'Saved {i} batches')
