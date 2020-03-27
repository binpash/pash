#!/bin/bash

# Identify top 10 terms in an input

N=10
IN=./input/1G.txt
OUT=./output/out.txt

cat $IN |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  sort |
  uniq -c |
  sort -rn |
  sed ${N}q > $OUT

