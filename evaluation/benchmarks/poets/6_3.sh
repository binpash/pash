#!/bin/bash
# tag: words_no_vowels
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/6_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c > ${OUT}/${input}.out
done

#for output in $(ls ${OUT} | sed "s;^;$OUT;")
#do
#    cat $output
#done

echo 'done';
rm -rf "${OUT}"
