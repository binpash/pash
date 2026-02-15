#!/bin/bash
#how the set-difference between two streams (i.e., elements in the first that are not in the second).
# https://stackoverflow.com/questions/2509533/bash-linux-set-difference-between-two-text-files

cd "$(dirname "$0")" || exit 1

# SIZE=500M
# IN="input/$SIZE.txt"

s1=${OUT}s1
s2=${OUT}s2

cat $IN |
    cut -d ' ' -f 1 |
    tr "[:lower:]" "[:upper:]" |
    sort > $s1

cat $IN |
    cut -d ' ' -f 1 |
    sort > $s2

comm -23 $s1 $s2 > ${OUT}stdout.txt

# rm s1 s2

