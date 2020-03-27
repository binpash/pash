#!/bin/bash

# Substring search  of  a complex  string, using  backtracking  (a superset  of
# Thompson's high-perf NFA used in GNU)

# Data: ./input/1G.txt

IN=./input/1G.txt
OUT=./output/out.txt

cat $IN |
  grep '[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}' > $OUT


