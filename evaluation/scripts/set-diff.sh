#!/bin/bash

# Show the set-difference between two streams (i.e., elements in the first that are not in the second).

# For default data, it uses the current set 

p="../../distributability/c_stats/"
A=${1:-${p}posix.txt}
B=${2:-${p}coreutils.txt}

IN=./input/i1G.txt
OUT=./output/out.txt

comm -23 <(cut -d ' ' -f 1 $IN | sort ) <( cut -d ' ' -f 1 $IN | sort) > $OUT
