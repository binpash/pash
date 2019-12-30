#!/bin/bash

# Find all 2-grams
# (old version, using intermediary files)

IN=./input/i1G.txt
OUT=./output/out.txt

cat $IN | 
  tr -cs A-Za-z '\n' |
  tr A-Z a-z > tokens.txt && 
tail +2 tokens.txt > next.txt &&
paste tokens.txt next.txt > bigrams.txt &&
cat bigrams.txt |
  sort | 
  uniq > $OUT
rm next.txt bigrams.txt
