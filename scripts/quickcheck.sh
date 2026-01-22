#!/bin/bash

# Get repo root from PASH_TOP (which points to src/pash)
REPO_ROOT="${PASH_TOP}/../.."
cd "$REPO_ROOT"

echo  confirms the necessary components for running the artifact
echo
echo Git commit ID: $(git rev-parse --short HEAD)
echo \$PASH_TOP: $(echo $PASH_TOP)
echo pash command: $(which pash 2>/dev/null || echo "not in PATH - activate virtual environment")

echo
pash --help

echo "Testing graph generation"
pash -c 'echo Pash Installation is complete!'
