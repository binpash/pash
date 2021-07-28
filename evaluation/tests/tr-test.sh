#!/bin/bash

## This test contains all occurences of tr (to test the annotation)

FILE="$PASH_TOP/evaluation/tests/input/1M.txt"

cat $FILE | tr -d ','
cat $FILE | tr '[A-Z]' '[a-z]'
cat $FILE | tr -s ' ' '\n'
cat $FILE | tr '[a-z]' 'P'
cat $FILE | tr -c "[a-z][A-Z]" '\n'
cat $FILE | tr ' ' '\n'
cat $FILE | tr '[a-z]' '\n'
## This is a bit tricky but `tr -d '\n'` is pure because after it is done there is only one line.
cat $FILE | tr -d '\n' | grep "the"
cat $FILE | tr -c '[A-Z]' '\n'
cat $FILE | tr " " " "
cat $FILE | tr -cs A-Za-z '\n'
cat $FILE | tr A-Z a-z
cat $FILE | tr -d '[:punct:]'
cat $FILE | tr [:lower] [:upper]
cat $FILE | tr [:lower:] [:upper:]
cat $FILE | tr -s ' '
cat $FILE | tr -s ' \n'
cat $FILE | tr -d '\012' | sort
