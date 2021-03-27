#!/bin/bash

# Identify top 10 terms in an input
# Edit: replaced `sed ${N}q` with `head -n$N`

N=10
IN=./input/1G.txt
OUT=./output/out.txt

cat $IN |
  sed 's/[^a-zA-Z0-9]/ /g' |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  sort |
  uniq -c |
  sort -rn |
  head -n$N > $OUT

