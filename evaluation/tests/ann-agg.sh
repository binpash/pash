#!/bin/bash

## Test contains command aliases with annotations that point to custom aggregators

FILE="../evaluation/tests/input/1M.txt"

test-one() {
  cat
}

test-two() {
  cat
}

cat $FILE | test-one | test-two
