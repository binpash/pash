from sklearn.linear_model import _logistic
import sys
from pickle import dumps, load

[solver, penalty, dual] = load(sys.stdin)
sys.stdout.buffer.write(dumps(_logistic._check_solver()))