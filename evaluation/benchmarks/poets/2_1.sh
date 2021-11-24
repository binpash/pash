#!/bin/bash
# tag: merge_upper
# set -e

# Merge upper and lower counts
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
# ls ${IN} | sed "s;^;$IN;"| xargs cat | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c

OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/2_1/}
mkdir -p "$OUT"

for input in $(ls ${IN} | sed "s;^;$IN;")
do
    cat "$input" | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c > "${OUT}/${input}"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

# Potentially rm -rf OUT
# rm 