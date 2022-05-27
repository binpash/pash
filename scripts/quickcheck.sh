#!/bin/bash

cd $PASH_TOP

echo  confirms the necessary components for running the artifact
echo
echo Git commit ID: $(git rev-parse --short HEAD)
echo \$PASH_TOP: $(echo $PASH_TOP)
echo pash executable: $PASH_TOP/pa.sh

echo
$PASH_TOP/pa.sh --help

echo "Testing graph generation"
$PASH_TOP/pa.sh -c 'echo Pash Installation is complete!'
