#!/bin/bash

## Running it with PaSh:
## time ./pa.sh -w 4 -d 1 --output_time evaluation/scripts/innefficient_auto_split.sh
##
## is slower than running it with bash:
## time ./evaluation/scripts/innefficient_auto_split.sh
##
## because the script doesn't do a lot of processing so 

FILE="$PASH_TOP/evaluation/scripts/input/1G.txt"
cat $FILE | sed 1d | grep 'Bell' | cut -f 2 | wc -l

## If instead we run the following, we get the expected results
# cat $FILE $FILE | grep 'Bell' | cut -f 2 | wc -l
