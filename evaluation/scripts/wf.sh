#!/bin/bash

# Word frequencies:

IN=./input/1G.txt
OUT=./output/out.txt

cat $IN |
  tr -cs A-Za-z'\n' |
  tr A-Z a-z |
  sort |
  uniq -c |
  sort > $OUT


