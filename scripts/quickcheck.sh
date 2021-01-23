#!/bin/bash

echo  confirms the necessary components for running the artifact
echo Git commit ID: $(git rev-parse --short HEAD)
echo \$PASH_TOP: $(echo $PASH_TOP)
echo pash executable: $(which pash)
echo Small input: $(wc -l ../evaluation/scripts/input/1M.txt)

../pa.sh --help


