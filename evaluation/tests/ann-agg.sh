#!/bin/bash

## Test contains command aliases with annotations that point to custom aggregators

FILE="$PASH_TOP/evaluation/tests/input/1M.txt"

test_one() {
  cat
}

test_two() {
  cat
}

cat $FILE | test_one | test_two
