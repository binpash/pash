from sklearn.linear_model import _logistic
import sys
import numpy as np
import pickle

model_file, y_file = sys.argv[1:]

with open(y_file, 'rb') as file:
    y = pickle.load(file)

try:
    # Can probably optimize, only opening the file once
    _logistic.check_classification_targets(y)
    with open(model_file, 'rb') as file:
        model = pickle.load(file)
        model.classes_ = np.unique(y)
    with open(model_file, 'wb') as file:
        pickle.dump(model, file)
except:
    exit(1)
