import sys
from torch import load, save

optimizer = load(sys.argv[1])
optimizer.zero_grad()
save(optimizer, 'optimizer.pt')