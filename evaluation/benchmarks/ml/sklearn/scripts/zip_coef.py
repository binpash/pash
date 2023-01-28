import sys
import numpy as np
import pickle

with open(sys.argv[1], 'r+b') as file:
    model = pickle.load(file)
    fold_coefs_ = pickle.load(sys.stdin.buffer)

    fold_coefs_, _, n_iter_ = zip(*fold_coefs_)

    sys.stdout.buffer.write(pickle.dumps(fold_coefs_))

    model.n_iter = np.asarray(n_iter_, dtype=np.int32)[:, 0]
    pickle.dump(model, file)
