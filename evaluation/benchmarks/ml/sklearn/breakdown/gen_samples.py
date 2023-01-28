from sklearn.model_selection import train_test_split
from sklearn.datasets import load_digits
import pickle

digits = load_digits()
data = train_test_split(digits.data, digits.target, test_size=0.2, random_state=0)

filenames = ['X_train', 'X_test', 'y_train', 'y_test']
files = [open(f'./tmp/{name}.obj', 'w+b') for name in filenames]
for datum, file in zip(data, files):
    pickle.dump(datum, file)
