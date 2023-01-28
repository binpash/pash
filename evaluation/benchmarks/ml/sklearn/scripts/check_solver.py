from sklearn.linear_model import _logistic
import sys
import pickle

with open(sys.argv[1], 'rb') as file:
    model = pickle.load(file)
    try:
        _logistic._check_solver(model.solver, model.penalty, model.dual)
        exit(0)
    except ValueError:
        exit(1)
