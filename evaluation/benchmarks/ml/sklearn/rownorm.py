from sklearn.linear_model import _logistic
import sys
from pickle import dumps, load

X = load(sys.stdin)
max_squared_sum = _logistic.row_norms(X, squared=True).max()
sys.stdout.buffer.write(dumps(max_squared_sum))
