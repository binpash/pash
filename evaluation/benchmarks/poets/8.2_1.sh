#!/bin/bash
# tag: vowel_sequences_gr_1K.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.2_1}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000"
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000" > "${OUT}/$(basename ${input})"
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
