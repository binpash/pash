#!/bin/bash
set -e

INTERMED_DIR=${INTERMED_DIR:-$PASH_TOP/evaluation/benchmarks/ml/intermed/}

TRAIN_LOADER=${TRAIN_LOADER:-$PASH_TOP/evaluation/benchmarks/ml/dataloader/train.pt}
PYTHON=${PYTHON:-`which python`}

$PYTHON network_generator.py $INTERMED_DIR
echo 'Finished generating network elements'

$PYTHON batchify.py $TRAIN_LOADER
echo 'Produced batches'