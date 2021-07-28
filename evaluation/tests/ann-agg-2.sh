#!/bin/bash

## Test contains command aliases with annotations that point to custom aggregators

FILE="$PASH_TOP/evaluation/tests/input/ab.txt"

test_uniq_1() {
  uniq
}

test_uniq_2() {
  uniq -c
}

cat $FILE | sort | test_uniq_1 | tr 'a' 'b' | test_uniq_2

