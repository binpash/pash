from sklearn.linear_model import _logistic
from pickle import load, dumps
import sys

reg = load(sys.stdin)
# X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, multi_class, 
# solver, penalty, sample_weight, n_threads
elem = []
for i in range(12):
    elem.append(load(sys.stdin))

path_func = _logistic.delayed(_logistic._logistic_regression_path)
fold_coefs = _logistic.Parallel(n_jobs=reg.n_jobs, verbose=reg.verbose, prefer=elem[5])(
    path_func(
        elem[0],
        elem[1],
        pos_class=class_,
        Cs=[elem[2]],
        l1_ratio=reg.l1_ratio,
        fit_intercept=reg.fit_intercept,
        tol=reg.tol,
        verbose=reg.verbose,
        solver=elem[8],
        multi_class=elem[7],
        max_iter=reg.max_iter,
        class_weight=reg.class_weight,
        check_input=False,
        random_state=reg.random_state,
        coef=warm_start_coef_,
        penalty=elem[9],
        max_squared_sum=elem[6],
        sample_weight=elem[10],
        n_threads=elem[11],
    )
    for class_, warm_start_coef_ in zip(elem[3], elem[4])
)
sys.stdout.buffer.write(dumps(fold_coefs))