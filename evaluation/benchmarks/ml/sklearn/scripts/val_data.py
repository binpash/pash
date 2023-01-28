import sys
import numpy as np
import pickle

with open(sys.argv[2], 'r+b') as X_file:
    X = pickle.load(X_file)
    with open(sys.argv[3], 'r+b') as y_file:
        y = pickle.load(y_file)
        with open(sys.argv[1], 'r+b') as model_file:
            model = pickle.load(model_file)
            try:
                X, y = model._validate_data(
                    X,
                    y,
                    accept_sparse="csr",
                    dtype=np.float64 if model.solver == 'lbfgs' else [np.float64, np.float32],
                    order="C",
                    accept_large_sparse=model.solver not in ["liblinear", "sag", "saga"],
                )

                pickle.dump(X, X_file)
                pickle.dump(y, y_file)
                pickle.dump(model, model_file)
                exit(0)
            except ValueError:
                exit(1)
