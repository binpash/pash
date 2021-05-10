#!/bin/bash

# Show the set-difference between two streams (i.e., elements in the first that are not in the second).

# For default data, it uses the current set 

# p="../../distributability/c_stats/"
# A=${1:-${p}posix.txt}
# B=${2:-${p}coreutils.txt}

# TODO: Find a good input set

mkfifo s1 s2

cat $IN |
    cut -d ' ' -f 1 |
    tr [:lower:] [:upper:] |
    sort > s1 &

cat $IN |
    cut -d ' ' -f 1 |
    sort > s2 &

comm -23 s1 s2

rm s1 s2
