PYTHON=${PYTHON:-`which python`}
DIR=$PASH_TOP/evaluation/benchmarks/ml/sklearn/breakdown

# takes solver, penalty, dual
# returns solver
$PYTHON $DIR/check_solver.py 
# takes model
# returns C_, penalty
$PYTHON $DIR/warning_checks.py 
# takes X, y, model, solver
# returns X, y
$PYTHON $DIR/val_data.py
# takes y 
# returns classes
$PYTHON $DIR/classes.py 
# takes model, solver, classes
# returns multi_class
$PYTHON $DIR/check_multiclass.py 
# takes X
# returns max_squared_sum
$PYTHON $DIR/rownorm.py 
# takes classes
# returns n_classes, classes reshaped
$PYTHON $DIR/reshape_classes.py
# takes model, multi_class, n_classes
# returns warm_start_coef
$PYTHON $DIR/warm_start.py 
# takes model, X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, 
# multi_class, solver, penalty, sample_weight, n_threads
# returns fold_coef
$PYTHON $DIR/fold_coef.py
# takes fold_coef
# returns fold_coef, n_iter
$PYTHON $DIR/zip_coef.py 
# takes X, multi_class, fold_coef, fit_intercept
# returns coef
$PYTHON $DIR/multiclass_coef.py 
# takes fit_intercept, coef 
# returns intercept and coef or zeroes
$PYTHON $DIR/fit_intercept.py 
# change model according to results
$PYTHON $DIR/fit_model.py 
