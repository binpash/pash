PASH_TOP=~/repos/research/pash # TODO: remove before commit

PYTHON=${PYTHON:-`which python`}
DIR=$PASH_TOP/evaluation/benchmarks/ml/sklearn/breakdown

# We can technically write to file via bash shovel, but for now, I write inside 
# the python script
MODEL=$DIR/tmp/model.obj
X=$DIR/tmp/X_train.obj
y=$DIR/tmp/y_train.obj
CLASSES=$DIR/tmp/classes.obj
DUAL=false # should be converted to bool inside script
MAX_SQ_SUM=$DIR/tmp/max_squared_sum.obj

# Generating model & samples
$PYTHON $DIR/gen_model.py
$PYTHON $DIR/gen_samples.py

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $DIR/check_solver.py $MODEL
$PYTHON $DIR/val_data.py $X $y $MODEL
$PYTHON $DIR/classes.py $y # This should return a classes with just the unique classes in y
MULTICLASS=$($PYTHON $DIR/check_multiclass.py $MODEL $CLASSES)

# Calculations functions
$PYTHON $DIR/rownorm.py $X > max_squared_sum.txt
$PYTHON $DIR/reshape_classes.py < classes.txt | # pipes n_classes, reshaped classes
$PYTHON $DIR/warm_start.py $MODEL | # pipes coefficients
$PYTHON $DIR/fold_coef.py $X $y (C_ from warning_checks.py) $PREFER $MULTICLASS $SOLVER $PENALTY < max_squared_sum.txt |
$PYTHON $DIR/zip_coef.py |
$PYTHON $DIR/multiclass_coef.py $X $MULTICLASS $FIT_INTERCEPT |
$PYTHON $DIR/fit_intercept.py $FIT_INTERCEPT |
$PYTHON $DIR/fit_model.py $MODEL n_iter.txt reshaped_classes.txt > result


# # takes solver, penalty, dual
# # returns solver
# $PYTHON $DIR/check_solver.py 

# # takes model
# # returns C_, penalty
# $PYTHON $DIR/warning_checks.py 

# # takes X, y, model, solver
# # returns X, y
# $PYTHON $DIR/val_data.py

# # takes y 
# # returns classes
# $PYTHON $DIR/classes.py 

# # takes X
# # returns max_squared_sum
# $PYTHON $DIR/rownorm.py

# # takes model, solver, classes
# # returns multi_class
# $PYTHON $DIR/check_multiclass.py 

# # TODO: fix implementation
# # takes classes
# # returns n_classes, classes reshaped
# $PYTHON $DIR/reshape_classes.py

# # takes model, multi_class, n_classes
# # returns warm_start_coef
# $PYTHON $DIR/warm_start.py 

# # takes model, X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, 
# # multi_class, solver, penalty, sample_weight, n_threads
# # returns fold_coef
# $PYTHON $DIR/fold_coef.py

# # takes fold_coef
# # returns fold_coef, n_iter
# $PYTHON $DIR/zip_coef.py 

# # takes X, multi_class, fold_coef, fit_intercept
# # returns coef
# $PYTHON $DIR/multiclass_coef.py 

# # takes fit_intercept, coef 
# # returns intercept and coef or zeroes
# $PYTHON $DIR/fit_intercept.py 

# # takes model, coef, n_iter, classes, intercept
# $PYTHON $DIR/fit_model.py 
