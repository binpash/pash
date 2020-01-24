#!/bin/bash

# Compares two streams element by element

# Data: ./input/input.txt

IN=./input/i1G.txt
OUT=./output/out.txt

# https://crashingdaily.wordpress.com/2008/03/06/diff-two-stdout-streams/
# This would work with coreutils.
shuf() { awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' "$@" | sort -k1,1n | cut -d ' ' -f2-; }
gen() { 
  cat $IN
#   head -n 2 
}
# alias gen='cat ./input/input.txt'


mkfifo s1 s2

cat $IN |
# head -n 2 | # small input
  shuf |
  sort |
  tr [:lower:] [:upper:] > s1 & 

cat $IN |
# head -n 2 | # small input
  shuf | 
  sort | 
  tr [:upper:] [:lower:] > s2 &

diff -B s1 s2 > $OUT
rm s1 s2
