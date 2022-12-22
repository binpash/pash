import sys
import numpy as np
from pickle import dumps, load

X = load(sys.stdin)
multi_class = load(sys.stdin)
fold_coefs_ = load(sys.stdin)
fit_intercept = load(sys.stdin)
n_classes = load(sys.stdin)


n_features = X.shape[1]
if multi_class == "multinomial":
    coef_ = fold_coefs_[0][0]
else:
    coef_ = np.asarray(fold_coefs_)
    coef_ = coef_.reshape(
        n_classes, n_features + int(fit_intercept)
    )

sys.stdout.buffer.write(dumps(coef_))