import sys
import numpy as np
import pickle

with open(sys.argv[1], 'r+b') as file1:
    model = pickle.load(file1)
    with open('./tmp/fold_coef.obj', 'r+b') as file2:
        fold_coefs_ = pickle.load(file2)

    fold_coefs_, _, n_iter_ = zip(*fold_coefs_)

    model.n_iter = np.asarray(n_iter_, dtype=np.int32)[:, 0]
    pickle.dump(model, file1)

with open('./tmp/fold_coef.obj', 'w+b') as file:
    pickle.dump(fold_coefs_, file)
