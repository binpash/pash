#!/bin/bash

# Word frequencies:

IN=./input/1G.txt
OUT=./output/out.txt

cat $IN |
  sed 's/[^a-zA-Z0-9]/ /g' |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  sort |
  uniq -c |
  sort > $OUT


