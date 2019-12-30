#!/bin/bash

# Count all words in a stream.

IN=./input/i1G.txt
OUT=./output/out.txt

cat $IN |
  wc -w > $OUT
