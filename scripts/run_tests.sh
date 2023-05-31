#!/usr/bin/env bash

set -x e

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

echo "Running intro tests..."
cd "$PASH_TOP/evaluation/intro"
./test.sh

echo "Running interface tests..."
cd "$PASH_TOP/evaluation/tests/interface_tests"
./run.sh

echo "Running compiler tests..."
cd "$PASH_TOP/evaluation/tests/"
./test_evaluation_scripts.sh
