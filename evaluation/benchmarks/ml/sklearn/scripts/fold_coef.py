from sklearn.linear_model import _logistic
import pickle
import sys

file_args = []
for file_name in sys.argv[1:8]:
    with open(file_name, 'rb') as file:
        file_args.append(pickle.load(file))
[model, X, y, C_, classes, warm_start_coef, max_squared_sum] = file_args

multi_class = sys.argv[8]
penalty = sys.argv[9]

path_func = _logistic.delayed(_logistic._logistic_regression_path)

# The SAG solver releases the GIL so it's more efficient to use
# threads for this solver.
if model.solver in ["sag", "saga"]:
    prefer = "threads"
else:
    prefer = "processes"

# TODO: Refactor this to avoid joblib parallelism entirely when doing binary
# and multinomial multiclass classification and use joblib only for the
# one-vs-rest multiclass case.
if (
    model.solver in ["lbfgs", "newton-cg", "newton-cholesky"]
    and len(classes) == 1
    and _logistic.effective_n_jobs(model.n_jobs) == 1
):
    # In the future, we would like n_threads = _openmp_effective_n_threads()
    # For the time being, we just do
    n_threads = 1
else:
    n_threads = 1

print(classes)

fold_coefs_ = _logistic.Parallel(n_jobs=model.n_jobs, verbose=model.verbose, prefer=prefer)(
    path_func(
        X,
        y,
        pos_class=class_,
        Cs=[C_],
        l1_ratio=model.l1_ratio,
        fit_intercept=model.fit_intercept,
        tol=model.tol,
        verbose=model.verbose,
        solver=model.solver,
        multi_class=multi_class,
        max_iter=model.max_iter,
        class_weight=model.class_weight,
        check_input=False,
        random_state=model.random_state,
        coef=warm_start_coef_,
        penalty=penalty,
        max_squared_sum=max_squared_sum,
        sample_weight=None, # Keep it as None for now.
        # For some reason, the n_threads breaks the whole process
        # n_threads=n_threads,
    )
    for class_, warm_start_coef_ in zip(classes, warm_start_coef)
)

with open('./tmp/fold_coef.obj', 'w+b') as file:
    pickle.dump(fold_coefs_, file)