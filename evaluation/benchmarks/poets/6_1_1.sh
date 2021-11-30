#!/bin/bash
# tag: uppercase_by_token
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/6_1_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^[A-Z]'

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat "$IN/$input" | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^[A-Z]' > "${OUT}/${input}.out"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
