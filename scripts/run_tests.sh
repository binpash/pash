#!/usr/bin/env bash

set -x e

# Source utils for helper functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
. "$SCRIPT_DIR/utils.sh"

# Always set PASH_TOP to the correct location (src/pash)
# This overrides any stale values from old PaSh installations
export PASH_TOP="$REPO_ROOT/src/pash"
# some tests expect TMPDIR to be set instead of using a default
export TMPDIR="${TMPDIR:-/tmp/}"

# dash only
echo "Running intro tests..."
cd "$REPO_ROOT/evaluation/intro"
./test.sh

# dash and bash
echo "Running interface tests..."
cd "$REPO_ROOT/evaluation/tests/interface_tests"
./run.sh

# dash and bash
echo "Running compiler tests..."
cd "$REPO_ROOT/evaluation/tests/"
./test_evaluation_scripts.sh

# bash only
# Note: These should only run nightly and not in the main workflow since they take a lot of time
# echo "Running bash interface tests..."
# cd "$REPO_ROOT/evaluation/tests/interface_tests"
# ./run_bash_tests.sh
