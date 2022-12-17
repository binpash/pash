from sklearn.linear_model import _logistic
from pickle import load, dumps
import sys

reg = load(sys.stdin)
(X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, multi_class, solver, penalty, sample_weight, n_threads) = load(sys.stdin)

path_func = _logistic.delayed(_logistic._logistic_regression_path)
fold_coefs = _logistic.Parallel(n_jobs=reg.n_jobs, verbose=reg.verbose, prefer=prefer)(
    path_func(
        X,
        y,
        pos_class=class_,
        Cs=[C_],
        l1_ratio=reg.l1_ratio,
        fit_intercept=reg.fit_intercept,
        tol=reg.tol,
        verbose=reg.verbose,
        solver=solver,
        multi_class=multi_class,
        max_iter=reg.max_iter,
        class_weight=reg.class_weight,
        check_input=False,
        random_state=reg.random_state,
        coef=warm_start_coef_,
        penalty=penalty,
        max_squared_sum=max_squared_sum,
        sample_weight=sample_weight,
        n_threads=n_threads,
    )
    for class_, warm_start_coef_ in zip(classes_, warm_start_coef)
)
sys.stdout.buffer.write(dumps(fold_coefs))