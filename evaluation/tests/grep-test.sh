#!/bin/bash

## This test contains all occurences of tr (to test the annotation)

FILE="$PASH_TOP/evaluation/tests/input/1M.txt"

cat $FILE | grep "the"
cat $FILE | grep -c "the"
