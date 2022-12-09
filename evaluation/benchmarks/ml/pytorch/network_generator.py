import sys
from torch import nn, optim, save

# Characteristics inherent to the dataset
INPUT_SIZE = 28**2 # Images are 28px * 28px
OUTPUT_SIZE  = 10

# Neural network characteristics
HIDDEN_LAYER_SIZES = [128, 64]
BATCH_SIZE = 32

FILENAMES = ['model.pt', 'criterion.pt', 'optimizer.pt']

def generate_network():
    # Set up the model
    model = nn.Sequential(nn.Linear(INPUT_SIZE, HIDDEN_LAYER_SIZES[0]), 
                        nn.ReLU(),
                        nn.Linear(HIDDEN_LAYER_SIZES[0], HIDDEN_LAYER_SIZES[1]),
                        nn.ReLU(),
                        nn.Linear(HIDDEN_LAYER_SIZES[1], OUTPUT_SIZE),
                        nn.LogSoftmax(dim=1))
    # negative log-likelihood loss
    criterion = nn.NLLLoss()
    optimizer = optim.SGD(model.parameters(), lr=0.004, momentum=0.9)

    return model, criterion, optimizer

def save_network(model, criterion, optimizer, dest_dir):
    for i, elem in enumerate([model, criterion, optimizer]):
        save(elem, dest_dir + FILENAMES[i])
        print(dest_dir + FILENAMES[i])

def main():
    if len(sys.argv) != 2:
        print('Incorrect Usage')
        exit(1)
    
    model, criterion, optimizer = generate_network()
    save_network(model, criterion, optimizer, dest_dir=sys.argv[1])

if (__name__ == '__main__'):
    main()
