from sklearn.linear_model import _logistic
import sys
from pickle import dumps, load

reg = load(sys.stdin)
solver = load(sys.stdin)
classes = load(sys.stdin)

multi_class = _logistic._check_multi_class(reg.multi_class, solver, len(reg.classes_))
sys.stdout.buffer.write(dumps(multi_class))