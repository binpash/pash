#!/bin/bash

# Run tests for the PaSh preprocessor
#
# Usage:
#   ./run_tests.sh           # Run all tests
#   ./run_tests.sh -v        # Run all tests with verbose output
#   ./run_tests.sh TestClass # Run specific test class
#   ./run_tests.sh TestClass.test_method  # Run specific test method

set -e

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASH_TOP="$(dirname "$SCRIPT_DIR")"

# Set up environment
export PASH_TOP
export PASH_TMP_PREFIX="${PASH_TMP_PREFIX:-$(mktemp -d /tmp/pash_test_XXXXXXX)/}"
export PASH_BASH_VERSION="${PASH_BASH_VERSION:-5 2 21}"
export PYTHONPATH="$SCRIPT_DIR:$PASH_TOP/compiler:$PYTHONPATH"

# Activate virtual environment if it exists
if [ -f "$PASH_TOP/python_pkgs/bin/activate" ]; then
    source "$PASH_TOP/python_pkgs/bin/activate"
fi

# Change to preprocessor directory
cd "$SCRIPT_DIR"

# Run tests
if [ $# -eq 0 ]; then
    # Run all tests
    python3 -m unittest discover -s shell_ast -p "test_*.py" -v
else
    # Run specified tests
    python3 -m unittest "shell_ast.test_walk_preprocess.$@"
fi

# Cleanup
if [ -d "$PASH_TMP_PREFIX" ]; then
    rm -rf "$PASH_TMP_PREFIX"
fi

echo ""
echo "All tests passed!"
