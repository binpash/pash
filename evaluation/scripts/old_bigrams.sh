#!/bin/bash

# Find all 2-grams
# (old version, using intermediary files)

cat ./input.txt | 
  tr -cs A-Za-z '\n' |
  tr A-Z a-z > tokens.txt && 
tail +2 tokens.txt > next.txt &&
paste tokens.txt next.txt > bigrams.txt &&
cat bigrams.txt |
sort | 
uniq > results
