from sklearn.linear_model import _logistic
import sys
import numpy as np
import pickle

with open(sys.argv[2], 'rb') as file:
    y = pickle.load(file)

try:
    _logistic.check_classification_targets(y)
    with open(sys.argv[1], 'rb') as file:
        model = pickle.load(file)
        model.classes_ = np.unique(y)
    with open(sys.argv[1], 'wb') as file:
        pickle.dump(model, file)
except:
    exit(1)
