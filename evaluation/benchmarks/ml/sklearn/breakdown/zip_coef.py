import sys
import numpy as np
from pickle import dumps, load

fold_coefs_ = load(sys.stdin)

fold_coefs_, _, n_iter_ = zip(*fold_coefs_)

sys.stdout.buffer.write(dumps(fold_coefs_))
sys.stdout.buffer.write(dumps(np.asarray(n_iter_, dtype=np.int32)[:, 0]))
