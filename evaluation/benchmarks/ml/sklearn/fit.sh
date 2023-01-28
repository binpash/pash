PASH_TOP=~/repos/research/pash # TODO: remove

PYTHON=${PYTHON:-`which python`}
DIR=$PASH_TOP/evaluation/benchmarks/ml/sklearn
TMP=$DIR/tmp
SCRIPTS=$DIR/scripts

# We can technically write to file via bash shovel, but for now, I write inside 
# the python script
MODEL=$TMP/model.obj
X=$TMP/X_train.obj
y=$TMP/y_train.obj
CLASSES=$TMP/classes.obj
DUAL=false # should be converted to bool inside script
MAX_SQ_SUM=$TMP/max_squared_sum.obj
WARM_COEF=$TMP/warm_start_coef.obj
C_=$TMP/C_.obj

# Generating model & samples
$PYTHON $SCRIPTS/gen_model.py
$PYTHON $SCRIPTS/gen_samples.py

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $SCRIPTS/check_solver.py $MODEL
penalty=$($PYTHON $SCRIPTS/penalty.py $MODEL)
$PYTHON $SCRIPTS/val_data.py $X $y $MODEL
$PYTHON $SCRIPTS/classes.py $y # This should return a classes with just the unique classes in y
multiclass=$($PYTHON $SCRIPTS/check_multiclass.py $MODEL $CLASSES)

# Calculations functions
$PYTHON $SCRIPTS/rownorm.py $X
n_classes=$($PYTHON $SCRIPTS/reshape_classes.py $CLASSES)
$PYTHON $SCRIPTS/warm_start.py $MODEL $multiclass $n_classes | # pipes coefficients
$PYTHON $SCRIPTS/fold_coef.py $MODEL $X $y $C_ $CLASSES $WARM_COEF $MAX_SQ_SUM $multiclass $penalty |
$PYTHON $SCRIPTS/zip_coef.py $MODEL |
$PYTHON $SCRIPTS/adjust_coef.py $MODEL $X $multiclass $n_classes


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