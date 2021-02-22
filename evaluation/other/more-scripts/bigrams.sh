#!/bin/bash

# Find all 2-grams in a piece of text

IN=./input/1G.txt
OUT=./output/out.txt

mkfifo s2
cat $IN |
  sed 's/[^a-zA-Z0-9]/ /g' |
# head -n 2 |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  tee s2 |
  tail +2 |
  paste s2 - |
  sort |
  uniq > $OUT
rm s2


