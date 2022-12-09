import torch
import sys

def train_batch(batch: torch.Tensor, labels: torch.Tensor) -> float:
    """ Trains the model for a batch of training data.

    Args:
        batch (torch.Tensor): the batch to train on
        labels (torch.Tensor): the label of the data

    Returns:
        float: the loss value from training this batch
    
    Currently relies on global variables for convenience. Should consider taking
    them as arguments if this were to be a separate script
    """
    # flattens size 28*28 matrix to length 784 array
    batch = batch.view(batch.shape[0], -1)

    # zero out gradients left over from previous iteration
    optimizer.zero_grad()
    # potential for pipelining
    output = model(batch)
    loss = criterion(output, labels)
    # back propagation
    loss.backward()
    # weight optimization
    optimizer.step()

    return loss.item()

def train_epoch(train_loader: torch.Tensor):
    running_loss = 0
    # iterates over each batch in the train dataset
    for images, labels in train_loader:
        running_loss += train_batch(batch=images, labels=labels)

def main():
    if len(sys.argv) != 5:
        print('Incorrect Usage')
        exit(1)
    global model, criterion, optimizer
    model, criterion, optimizer, train_loader = \
        [torch.load(path) for path in sys.argv[1:]]

    train_epoch(train_loader)
    
    torch.save(model, sys.argv[1])
    torch.save(criterion, sys.argv[2])
    torch.save(optimizer, sys.argv[3])
