#!/bin/bash
# https://crashingdaily.wordpress.com/2008/03/06/diff-two-stdout-streams/
# This would work with coreutils.
shuf() { awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' "$@" | sort -k1,1n | cut -d ' ' -f2-; }

alias gen='head -n 2 ./input/input.txt'
# alias gen='cat ./input/input.txt'

mkfifo s1 s2
gen | shuf | sort | tr [:lower:] [:upper:] > s1 & 
gen | shuf | sort | tr [:lower:] [:upper:] > s2 &
diff -B s1 s2
rm s1 s2
