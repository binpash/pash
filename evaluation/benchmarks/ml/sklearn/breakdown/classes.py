from sklearn.linear_model import _logistic
import sys
import numpy as np
from pickle import dumps, load

y = load(sys.stdin)

_logistic.check_classification_targets(y)
sys.stdout.buffer.write(dumps(np.unique(y)))