from sklearn.model_selection import train_test_split
from sklearn import datasets
import pickle
import sys
import pandas as pd

data_name = sys.argv[1]

if data_name == 'digit':
    raw_data = datasets.load_digits()
elif data_name == 'cover':
    raw_data = datasets.fetch_covtype()
    
data = train_test_split(raw_data.data, 
                        raw_data.target, 
                        test_size=0.2, 
                        random_state=0)
filenames = ['X_train', 'X_test', 'y_train', 'y_test']
files = [open(f'./tmp/{name}.obj', 'w+b') for name in filenames]
for datum, file in zip(data, files):
    pickle.dump(datum, file)
