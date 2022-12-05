import torch
import torchvision
import pathlib

from time import time
from torch import nn, optim
from torch.utils.data import DataLoader
from torchvision import transforms, datasets

# Characteristics inherent to the dataset
INPUT_SIZE = 28**2 # Images are 28px * 28px
OUTPUT_SIZE  = 10

# Neural network characteristics
HIDDEN_LAYER_SIZES = [128, 64]
BATCH_SIZE = 32
EPOCHS = 15

# Paths for storing datasets (to not re-download it every time)
DATASETS_PATH = pathlib.Path(__file__).parent.resolve() / 'datasets'
TRAINSET_PATH = DATASETS_PATH / 'train.txt' 
TESTSET_PATH = DATASETS_PATH / 'train.txt'
USE_GPU = True

def set_device() -> torch.device:
    """ Sets the device to train the model on between GPU and CPU

    Returns:
        torch.device: The appropriate device
    """
    if torch.cuda.is_available() and USE_GPU:
        print('Using Cuda')
        return torch.device('cuda:0')
    else:
        print('Using CPU')
        return torch.device('cpu')

def load_data(batch_size: int) -> tuple[DataLoader, DataLoader]:
    """downloads and processes the MNIST handwriting dataset

    Args:
        batch_size (int): the batch size to train with

    Returns:
        tuple[DataLoader, DataLoader]: a train dataset and a test dataset
    """
    transform = transforms.Compose([transforms.ToTensor(), 
                                    transforms.Normalize((0.5,), (0.5,))])

    # TODO: Generalize this for use on datasets not provided by torch
    train_set = datasets.MNIST(TRAINSET_PATH, download=True, train=True, 
                               transform=transform)
    test_set = datasets.MNIST(TESTSET_PATH, download=True, train=False, 
                              transform=transform)
    # Wrapper for DataSets for easy access to individual batches, etc
    train_loader = DataLoader(train_set, batch_size, shuffle=True)
    test_loader = DataLoader(test_set, batch_size, shuffle=True)

    return train_loader, test_loader

def train_batch(batch: torch.Tensor, labels: torch.Tensor, print_time) -> float:
    """ Trains the model for a batch of training data.

    Args:
        batch (torch.Tensor): the batch to train on
        labels (torch.Tensor): the label of the data

    Returns:
        float: the loss value from training this batch
    
    Currently relies on global variables for convenience. Should consider taking
    them as arguments if this were to be a separate script
    """
    batch_start_time = time()

    # flattens size 28*28 matrix to length 784 array
    batch = batch.view(batch.shape[0], -1)
    flatten_batch_time = time()

    # zero out gradients left over from previous iteration
    optimizer.zero_grad()
    zero_grad_time = time()

    # potential for pipelining
    output = model(batch)
    output_time = time()
    loss = criterion(output, labels)
    loss_time = time()

    # back propagation
    loss.backward()
    back_prop_time = time()

    # weight optimization
    optimizer.step()
    optimizer_time = time()

    if print_time:
        print(f'flatten_batch_time: {flatten_batch_time - batch_start_time}')
        print(f'zero_grad_time: {zero_grad_time - flatten_batch_time}')
        print(f'output_time: {output_time - zero_grad_time}')
        print(f'loss_time: {loss_time - output_time}')
        print(f'back_prop_time: {back_prop_time - loss_time}')
        print(f'optimizer_time: {optimizer_time - back_prop_time}')

    return loss.item()

def prepare_environment():
    # Set up the model
    device = set_device()
    model = nn.Sequential(nn.Linear(INPUT_SIZE, HIDDEN_LAYER_SIZES[0]), 
                        nn.ReLU(),
                        nn.Linear(HIDDEN_LAYER_SIZES[0], HIDDEN_LAYER_SIZES[1]),
                        nn.ReLU(),
                        nn.Linear(HIDDEN_LAYER_SIZES[1], OUTPUT_SIZE),
                        nn.LogSoftmax(dim=1)).to(device=device)
    # negative log-likelihood loss
    criterion = nn.NLLLoss()
    optimizer = optim.SGD(model.parameters(), lr=0.004, momentum=0.9)

    return device, model, criterion, optimizer

def train_epoch(n: int):
    running_loss = 0
    # iterates over each batch in the train dataset
    for images, labels in train_loader:
        running_loss += train_batch(batch=images.to(device=device), 
                                    labels=labels.to(device=device),
                                    print_time=False)
    print(f'Epoch {n} Loss: {running_loss / len(train_loader)}')

def calc_accuracy(model: nn.Sequential, test_loader: DataLoader, device: torch.device):
    correct_count, all_count = 0, 0
    for images, labels in test_loader:
        images = images.to(device=device)
        labels = labels.to(device=device)

        for i, image in enumerate(images):
            image = image.view(1, 784)
            
            # Don't want to change the gradient of model - just testing
            with torch.no_grad():
                log_probabilities = model(image)
            
            probabilities = torch.exp(log_probabilities)
            probability = list(probabilities.cpu().numpy()[0])
            predicted_label = probability.index(max(probability))
            true_label = labels.cpu().numpy()[i]
            if (true_label == predicted_label):
                correct_count += 1
            all_count += 1

    print(f'No. images tested: {all_count}')
    print(f'Accuracy: {correct_count / all_count}')

def main():
    device, model, criterion, optimizer = prepare_environment()
    train_loader, test_loader = load_data(BATCH_SIZE)
    time_0 = time()

    print('Begin Training...')
    for n in range(EPOCHS):
        train_epoch(n)
    print(f'Training time: {time() - time_0}')

    calc_accuracy(model, test_loader, device)