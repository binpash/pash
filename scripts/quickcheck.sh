#!/bin/bash

cd $PASH_TOP

echo  confirms the necessary components for running the artifact
echo
echo Git commit ID: $(git rev-parse --short HEAD)
echo \$PASH_TOP: $(echo $PASH_TOP)
echo pash executable: $PASH_TOP/pa.sh

# SMALL=$PASH_TOP/evaluation/scripts/input/3M.txt
echo 'Inputs:'
ls -lh $PASH_TOP/evaluation/scripts/input/* | awk '{print $5,$9}'


echo
$PASH_TOP/pa.sh --help
