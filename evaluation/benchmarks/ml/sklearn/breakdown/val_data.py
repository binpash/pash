import sys
from pickle import dumps, load
import numpy as np

X = load(sys.stdin)
y = load(sys.stdin)
reg = load(sys.stdin)
solver = load(sys.stdin)

X, y = reg._validate_data(
    X,
    y,
    accept_sparse="csr",
    dtype=np.float64 if solver == 'lbfgs' else [np.float64, np.float32],
    order="C",
    accept_large_sparse=solver not in ["liblinear", "sag", "saga"],
)