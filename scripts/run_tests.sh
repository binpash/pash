#!/usr/bin/env bash

set -x e

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# some tests expect TMPDIR to be set instead of using a default
export TMPDIR="${TMPDIR:-/tmp/}"

echo "Running intro tests..."
cd "$PASH_TOP/evaluation/intro"
./test.sh

echo "Running interface tests..."
cd "$PASH_TOP/evaluation/tests/interface_tests"
./run.sh

echo "Running compiler tests..."
cd "$PASH_TOP/evaluation/tests/"
./test_evaluation_scripts.sh
