import sys
from torch import load, save

optimizer = load(sys.argv[1])
optimizer.step()
save(optimizer, 'optimizer.pt')
