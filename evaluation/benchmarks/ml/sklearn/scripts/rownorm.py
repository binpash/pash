from sklearn.linear_model import _logistic
import sys
import pickle

with open(sys.argv[1], 'rb') as file:
    X = pickle.load(file)
    
max_squared_sum = _logistic.row_norms(X, squared=True).max()

with open('./tmp/max_squared_sum.obj', 'w+b') as file:
    pickle.dump(max_squared_sum, file)
