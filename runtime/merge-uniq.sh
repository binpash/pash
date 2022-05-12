#!/usr/bin/env bash

# This is how to merge results of `uniq -c`, contained in {1,2,3}.txt
# I am using 3 inputs to stress it works with more than just pairs:-)

A=${1:-1.txt}
B=${1:-2.txt}
C=${1:-3.txt}
awk '{ count[$2] += $1 } END { for(e in count) print count[e], e }' "$A" "$B" "$C"
