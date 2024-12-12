#!/usr/bin/env bash

set -x e

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

echo "Running intro tests..."
cd "$PASH_TOP/evaluation/intro"
./test.sh bash

echo "Running interface tests..."
cd "$PASH_TOP/evaluation/tests/interface_tests"
./run.sh bash

echo "Running compiler tests..."
cd "$PASH_TOP/evaluation/tests/"
./test_evaluation_scripts.sh bash

echo "Running bash interface tests..."
cd "$PASH_TOP/evaluation/tests/interface_tests"
./run_bash_tests.sh bash
