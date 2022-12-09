import sys
from torch import load, save

output = load(sys.argv[1])
labels = load(sys.argv[2])
criterion = load(sys.argv[3])
loss_path = sys.argv[4]

loss = criterion(output, labels)
loss.backward()
save(loss, loss_path)
