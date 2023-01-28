from sklearn.linear_model import _logistic
import sys
import pickle

with open(sys.argv[1], 'rb') as file:
    model = pickle.load(file)
with open(sys.argv[2], 'rb') as file:
    classes = pickle.load(file)

multi_class = _logistic._check_multi_class(model.multi_class, model.solver, len(classes))
print(multi_class)