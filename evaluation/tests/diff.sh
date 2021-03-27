#!/bin/bash

# Compares two streams element by element

# https://crashingdaily.wordpress.com/2008/03/06/diff-two-stdout-streams/
# This would work with coreutils.
# shuf() { awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' "$@" | sort -k1,1n | cut -d ' ' -f2-; }

mkfifo s1 s2

cat $IN |
  # shuf |
  tr [:lower:] [:upper:] |
  sort > s1 &

cat $IN |
  # shuf |
  tr [:upper:] [:lower:] |
  sort > s2 &

diff -B s1 s2
rm s1 s2
