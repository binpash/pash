#!/bin/bash

# Find all 2-grams in a piece of text
# FIXME: does not calculate frequencies

IN=./input/i10M.txt
OUT=./output/out.txt

## Currently our translation pass doesn't handle Defun

# mkfifo s2

# bigrams_aux()
# {
#     ## I am not sure if this reads the stdin of the function
#     tee s2 |
#         tail +2 |
#         paste s2 -
# }

cat $IN $IN |
# head -n 2 |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq > $OUT

# rm s2


