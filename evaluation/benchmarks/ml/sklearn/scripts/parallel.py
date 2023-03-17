from sklearn.linear_model import _logistic
import pickle
import sys

file_args = []
for file_name in sys.argv[1:7]:
    with open(file_name, 'rb') as file:
        file_args.append(pickle.load(file))
model, X, y, C_, warm_start_coef_, max_squared_sum = file_args

multi_class = sys.argv[7]
penalty = sys.argv[8]
class_ = int(sys.argv[9])

result = _logistic._logistic_regression_path(
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
    coef=None, #Leave as None for now, consider changing to warm_start_coef_ later
    penalty=penalty,
    max_squared_sum=max_squared_sum,
    sample_weight=None,
)

with open(f'./tmp/result_{class_}.obj', 'wb') as file:
    pickle.dump(result, file)