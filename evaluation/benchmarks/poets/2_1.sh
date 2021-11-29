#!/bin/bash
# tag: merge_upper
# set -e

# Merge upper and lower counts
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/2_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
# ls ${IN} | sed "s;^;$IN;"| xargs cat | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c

for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c > "${OUT}/$(basename ${input})"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
