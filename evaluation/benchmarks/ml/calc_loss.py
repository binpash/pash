import sys
from torch import load, save

output = load(sys.argv[1])
labels = load(sys.argv[2])
criterion = load(sys.argv[3])

save(criterion(output, labels).backward(), 'loss.pt')
