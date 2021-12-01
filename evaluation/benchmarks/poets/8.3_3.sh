#!/bin/bash
# tag: compare_exodus_genesis.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/benchmarks/poets/input/exodus}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u >  ${OUT}1.types
#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > ${OUT}2.types
## FIXME do we really need the same thing twice?
#sort ${OUT}1.types ${OUT}2.types ${OUT}2.types | uniq -c | head
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.2_1}
mkdir -p "$OUT"
ENTRIES=${ENTRIES:-1000}
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}1.types
    tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > ${OUT}/${input}2.types
    sort ${OUT}/${input}1.types ${OUT}/${input}2.types ${OUT}/${input}2.types | uniq -c | head
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
