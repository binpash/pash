import sys
import torch

from pathlib import Path
from torch.utils.data import DataLoader
from torchvision import transforms, datasets

# Paths for storing datasets (to not re-download it every time)
DATASETS_PATH = Path(__file__).parent.resolve() / 'datasets'
TRAINSET_PATH = DATASETS_PATH / 'train.txt' 
TESTSET_PATH = DATASETS_PATH / 'train.txt'

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
    # TODO: Fix shuffle seed to ensure same dataloader is generated
    train_loader = DataLoader(train_set, batch_size, shuffle=True)
    test_loader = DataLoader(test_set, batch_size, shuffle=True)

    return train_loader, test_loader

def main():
    if len(sys.argv) != 3:
        print('Incorrect Usage')
        print('Correct Use: data.py <batch-size> <destination-dir>')
        exit(1)

    # TODO: gracefully exit with try-catch

    batch_size = int(sys.argv[1])
    dest = sys.argv[2]
    train_dest = dest + 'train.pt'
    test_dest = dest + 'test.pt'

    train, test = load_data(batch_size)
    torch.save(train, train_dest)
    torch.save(test, test_dest)

    print(train_dest)
    print(test_dest)


if (__name__ == '__main__'):
    main()