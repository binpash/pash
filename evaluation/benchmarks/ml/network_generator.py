from torch import nn, optim

# Characteristics inherent to the dataset
INPUT_SIZE = 28**2 # Images are 28px * 28px
OUTPUT_SIZE  = 10

# Neural network characteristics
HIDDEN_LAYER_SIZES = [128, 64]
BATCH_SIZE = 32
EPOCHS = 15

def prepare_environment():
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