#!/bin/bash
# Find all 2-grams in a piece of text

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

. bi-gram.aux.sh

cat $IN |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq


