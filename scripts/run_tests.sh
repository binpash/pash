#!/usr/bin/env bash

set -x e

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# some tests expect TMPDIR to be set instead of using a default
export TMPDIR="${TMPDIR:-/tmp/}"

# dash only
echo "Running intro tests..."
cd "$PASH_TOP/evaluation/intro"
./test.sh

# dash and bash
echo "Running interface tests..."
cd "$PASH_TOP/evaluation/tests/interface_tests"
./run.sh

# dash and bash
echo "Running compiler tests..."
cd "$PASH_TOP/evaluation/tests/"
./test_evaluation_scripts.sh

# bash only
# Note: These should only run nightly and not in the main workflow since they take a lot of time
# echo "Running bash interface tests..."
# cd "$PASH_TOP/evaluation/tests/interface_tests"
# ./run_bash_tests.sh
