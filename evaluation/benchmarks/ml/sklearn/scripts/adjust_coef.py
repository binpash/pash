import sys
import pickle
import numpy as np

fold_coefs_ = pickle.load(sys.stdin.buffer)

multi_class = sys.argv[3]
n_classes = sys.argv[4]

with open(sys.argv[2], 'rb') as file:
    X = pickle.load(file)
    n_features = X.shape[1]

with open(sys.argv[1], 'rb') as file:
    model = pickle.load(file)
    if multi_class == "multinomial":
        model.coef_ = fold_coefs_[0][0]
    else:
        model.coef_ = np.asarray(fold_coefs_)
        model.coef_ = model.coef_.reshape(
            n_classes, n_features + int(model.fit_intercept)
        )
    
    if model.fit_intercept:
        model.intercept_ = model.coef_[:, -1]
        model.coef_ = model.coef_[:, :-1]
    else:
        model.intercept_ = np.zeros(n_classes)

with open('./trained_model.obj', 'w+b') as file:
    pickle.dump(model, file)