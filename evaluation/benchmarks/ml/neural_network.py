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
BATCH_SIZE = 256
EPOCHS = 15

# Paths for storing datasets (to not re-download it every time)
DATASETS_PATH = pathlib.Path(__file__).parent.resolve() / 'datasets'
TRAINSET_PATH = DATASETS_PATH / 'train.txt' 
TESTSET_PATH = DATASETS_PATH / 'train.txt'

USE_GPU = True

def set_device():
    if torch.cuda.is_available() and USE_GPU:
        print('Using Cuda')
        return torch.device('cuda:0')
    else:
        print('Using CPU')
        return torch.device('cpu')

DEVICE = set_device()

def load_data(batch_size: int):
    """downloads and processes the MNIST handwriting dataset

    Args:
        batch_size (int): the batch size to train with

    Returns:
        DataLoader: a train dataset and a test dataset
    """
    transform = transforms.Compose([transforms.ToTensor(), 
                                    transforms.Normalize((0.5,), (0.5,))])

    # TODO: Generalize this for use on datasets not provided by torch
    train_set = datasets.MNIST(TRAINSET_PATH, download=True, train=True, 
                               transform=transform)
    test_set = datasets.MNIST(TESTSET_PATH, download=True, train=False, 
                              transform=transform)

    train_loader = DataLoader(train_set, batch_size, shuffle=True)
    test_loader = DataLoader(test_set, batch_size, shuffle=True)

    return train_loader, test_loader

def train_batch(batch, labels):
    """ Trains the model for a batch of training data.

    Args:
        batch (_type_): the batch to train on
        labels: the label of the data

    Returns:
        float: the loss value from training this batch
    
    Currently relies on global variables for convenience. Should consider taking
    them as arguments if this were to be a separate script
    """

    # flattens size 28*28 matrix to length 784 array
    batch = batch.view(batch.shape[0], -1)

    # zero out gradients left over from previous iteration
    optimizer.zero_grad()

    output = model(batch)
    loss = criterion(output, labels)

    # back propagation
    loss.backward()
    # weight optimization
    optimizer.step()

    return loss.item()

# Set up the model
model = nn.Sequential(nn.Linear(INPUT_SIZE, HIDDEN_LAYER_SIZES[0]), 
                      nn.ReLU(),
                      nn.Linear(HIDDEN_LAYER_SIZES[0], HIDDEN_LAYER_SIZES[1]),
                      nn.ReLU(),
                      nn.Linear(HIDDEN_LAYER_SIZES[1], OUTPUT_SIZE),
                      nn.LogSoftmax(dim=1)).to(device=DEVICE)
# negative log-likelihood loss
criterion = nn.NLLLoss()
optimizer = optim.SGD(model.parameters(), lr=0.004, momentum=0.9)

train_loader, test_loader = load_data(BATCH_SIZE)

time_0 = time()

for n in range(EPOCHS):
    running_loss = 0
    # iterates over each batch in the train dataset
    for images, labels in train_loader:
        running_loss += train_batch(batch=images.to(device=DEVICE), 
                                    labels=labels.to(device=DEVICE))
    print(f'Epoch {n} Loss: {running_loss / len(train_loader)}')
print(f'Training time: {time() - time_0}')

correct_count, all_count = 0, 0
for images, labels in test_loader:
    images = images.to(device=DEVICE)
    labels = labels.to(device=DEVICE)

    for i, image in enumerate(images):
        image = image.view(1, 784)
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
