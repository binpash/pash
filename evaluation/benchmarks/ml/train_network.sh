#!/bin/bash
set -e

N_EPOCH=${N_EPOCH:-15}
INTERMED_DIR=${INTERMED_DIR:-$PASH_TOP/evaluation/benchmarks/ml/intermed/}
SCRIPTS_DIR=${SCRIPTS_DIR:-$PASH_TOP/evaluation/benchmarks/ml/python_pipeline/}

PYTHON=${PYTHON:-`which python`}

OPTIMIZER=${INTERMED_DIR}optimizer.pt
MODEL=${INTERMED_DIR}model.pt
CRITERION=${INTERMED_DIR}criterion.pt
OUTPUT=${INTERMED_DIR}output.pt
LOSS=${INTERMED_DIR}loss.pt
BATCHES_DIR=${INTERMED_DIR}batches/
LABELS_DIR=${INTERMED_DIR}labels/


N_BATCHES=`ls -v $BATCHES_DIR | tail -1 | tr -dc [:digit:]`

train_batch() {
    $PYTHON ${SCRIPTS_DIR}zero_grad.py $OPTIMIZER
    $PYTHON ${SCRIPTS_DIR}feed_batch_to_model.py $MODEL $1 $OUTPUT
    $PYTHON ${SCRIPTS_DIR}calc_loss.py $OUTPUT $2 $CRITERION $LOSS
    $PYTHON ${SCRIPTS_DIR}step_optimizer.py $OPTIMIZER
}

for epoch in $(seq 1 $N_EPOCH)
do   
    echo Epoch $epoch
    for i in $(seq 1 $N_BATCHES)
    do
        train_batch ${BATCHES_DIR}batch_$i.pt ${LABELS_DIR}labels_$i.pt
    done
done

echo 'done'