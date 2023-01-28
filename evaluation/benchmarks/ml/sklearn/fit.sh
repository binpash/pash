PASH_TOP=~/repos/research/pash # TODO: remove

PYTHON=${PYTHON:-`which python`}
DIR=$PASH_TOP/evaluation/benchmarks/ml/sklearn
TMP=$DIR/tmp
SCRIPTS=$DIR/scripts

# Ideally, we'll move on to piping rather than writing to a file
MODEL=$TMP/model.obj
X=$TMP/X_train.obj
y=$TMP/y_train.obj
CLASSES=$TMP/classes.obj
DUAL=false # should be converted to bool inside script
MAX_SQ_SUM=$TMP/max_squared_sum.obj
WARM_COEF=$TMP/warm_start_coef.obj
C_=$TMP/C_.obj

echo Destination: $1

# Generating model & samples
$PYTHON $SCRIPTS/gen_model.py
$PYTHON $SCRIPTS/gen_samples.py

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $SCRIPTS/check_solver.py $MODEL
penalty=$($PYTHON $SCRIPTS/penalty.py $MODEL)
$PYTHON $SCRIPTS/val_data.py $MODEL $X $y 
$PYTHON $SCRIPTS/classes.py $MODEL $y # This should return a classes with just the unique classes in y
multiclass=$($PYTHON $SCRIPTS/check_multiclass.py $MODEL)

# Calculations functions
$PYTHON $SCRIPTS/rownorm.py $X
n_classes=$($PYTHON $SCRIPTS/reshape_classes.py $MODEL $CLASSES)
$PYTHON $SCRIPTS/warm_start.py $MODEL $multiclass $n_classes | # pipes coefficients
$PYTHON $SCRIPTS/fold_coef.py $MODEL $X $y $C_ $CLASSES $WARM_COEF $MAX_SQ_SUM $multiclass $penalty |
$PYTHON $SCRIPTS/zip_coef.py $MODEL |
$PYTHON $SCRIPTS/adjust_coef.py $MODEL $X $multiclass $n_classes $1