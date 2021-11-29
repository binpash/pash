#!/bin/bash
# tag: sort_words_by_folding
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | sort -f 
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/3_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | sort -f > "${OUT}/$(basename ${input})"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
