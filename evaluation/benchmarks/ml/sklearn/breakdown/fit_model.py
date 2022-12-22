import sys
from pickle import dumps, load

reg = load(sys.stdin)
reg.coef_ = load(sys.stdin)
reg.n_iter_ = load(sys.stdin)
reg.classes_ = load(sys.stdin)
reg.intercept = load(sys.stdin)

if reg.fit_intercept:
    reg.coef_ = load(sys.stdin)

sys.stdout.buffer.write(dumps(reg))