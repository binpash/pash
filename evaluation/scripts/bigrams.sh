#!/bin/bash

# Find all 2-grams in a piece of text
# FIXME: does not calculate frequencies

IN=./input/1G.txt
OUT=./output/out.txt

mkfifo s2
cat $IN |
# head -n 2 |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  tee s2 |
  tail +2 |
  paste s2 - |
  sort |
  uniq > $OUT
rm s2


