#!/bin/bash

# Word frequencies:

IN=./input/i1G.txt
OUT=./output/out.txt

cat $IN |
  tr -cs A-Za-z'\n' |
  tr A-Z a-z |
  sort |
  uniq -c |
  sort > $OUT


