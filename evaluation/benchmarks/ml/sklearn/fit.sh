# TODO: figure out a way to set PASH_TOP in the python subprocess module
# so that I can delete the line below
PASH_TOP=~/repos/research/pash 

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

# TODO: Try this out on a larger dataset
# TODO: Benchmark each phase

# Generating model & samples
$PYTHON $SCRIPTS/gen_model.py $3
$PYTHON $SCRIPTS/gen_samples.py $2

echo Step 1 end

# Validity checking functions
# These functions just check to make sure that the input is valid. 
# If not they will raise an error. Otherwise, they do not mutate the data.
$PYTHON $SCRIPTS/check_solver.py $MODEL
penalty=$($PYTHON $SCRIPTS/penalty.py $MODEL)
$PYTHON $SCRIPTS/val_data.py $MODEL $X $y 
$PYTHON $SCRIPTS/classes.py $MODEL $y # This should return a classes with just the unique classes in y
multiclass=$($PYTHON $SCRIPTS/check_multiclass.py $MODEL)

echo Step 2 end

# TODO: Benchmark each step of the pipeline
# Make a modified pipeline where each step writes its output to a file

# Calculations functions
$PYTHON $SCRIPTS/rownorm.py $X
echo rownorm
n_classes=$($PYTHON $SCRIPTS/reshape_classes.py $MODEL $CLASSES)
echo reshape
$PYTHON $SCRIPTS/warm_start.py $MODEL $multiclass $n_classes # pipes coefficients
echo warm_start

# Covtype dataset has 7 classes
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 1
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 2
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 3
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 4
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 5
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 6
$PYTHON $SCRIPTS/parallel.py $MODEL $X $y $C_ $WARM_COEF $MAX_SQ_SUM $multiclass $penalty 7
echo fold_coef

$PYTHON $SCRIPTS/zip_coef.py $MODEL
echo zip_coef
$PYTHON $SCRIPTS/adjust_coef.py $MODEL $X $multiclass $n_classes $1
echo adjust_coef

echo Step 3 end
