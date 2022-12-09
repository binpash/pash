import torch
from torch import nn
from torch.utils.data import DataLoader

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
