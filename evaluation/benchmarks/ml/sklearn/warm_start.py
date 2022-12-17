import sys
import numpy as np
from pickle import dumps, load

reg = load(sys.stdin)
multi_class = load(sys.stdin)
n_classes = load(sys.stdin)

if reg.warm_start:
    warm_start_coef = getattr(reg, "coef_", None)
else:
    warm_start_coef = None
if warm_start_coef is not None and reg.fit_intercept:
    warm_start_coef = np.append(
        warm_start_coef, reg.intercept_[:, np.newaxis], axis=1
    )

if multi_class == "multinomial":
    classes_ = [None]
    warm_start_coef = [warm_start_coef]
if warm_start_coef is None:
    warm_start_coef = [None] * n_classes

sys.stdout.buffer.write(dumps(warm_start_coef))