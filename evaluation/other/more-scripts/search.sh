#!/bin/bash

# Complicated grep expression

IN=./input/1G.txt # Change G to M for small input
OUT=./output/out.txt

cat $IN | 
  grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > $OUT
