#!/bin/bash
set -e

INTERMED_DIR=${INTERMED_DIR:-$PASH_TOP/evaluation/benchmarks/ml/pytorch/intermed/}

TRAIN_LOADER=${TRAIN_LOADER:-$PASH_TOP/evaluation/benchmarks/ml/pytorch/dataloader/train.pt}
PYTHON=${PYTHON:-`which python`}

$PYTHON network_generator.py $INTERMED_DIR
echo "Done generating network elements"

$PYTHON batchify.py $TRAIN_LOADER
echo "Done producing batches"