PYTHON=${PYTHON:-`which python`}

$PYTHON check_solver.py # returns solver
$PYTHON warning_checks # returns C_, penalty
$PYTHON val_data # returns X, y
$PYTHON classes.py # returns classes
$PYTHON check_multiclass.py # takes model, solver, classes; returns multi_class

$PYTHON rownorm.py # takes X; returns max_squared_sum
$PYTHON reshape_classes.py # takes classes; returns n_classes, classes reshaped
$PYTHON warm_start.py # takes model, multi_class, n_classes; returns warm_start_coef
# takes model, 
# (X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, 
# multi_class, solver, penalty, sample_weight, n_threads)
# returns fold_coef
$PYTHON fold_coef.py