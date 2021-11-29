#!/bin/bash 
#tag: count_trigrams.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/4_3b/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
#tail +2 ${OUT}.words > ${OUT}.nextwords
#tail +3 ${OUT}.words > ${OUT}.nextwords2
#paste ${OUT}.words ${OUT}.nextwords ${OUT}.nextwords2 |
#sort | uniq -c  > ${OUT}.trigrams
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr -sc '[A-Z][a-z]' '[\012*]' > "${OUT}/$(basename ${input}).words"
    tail +2 "${OUT}/$(basename ${input}).words" > "${OUT}/$(basename ${input}).nextwords"
    tail +2 "${OUT}/$(basename ${input}).words" > "${OUT}/$(basename ${input}).nextwords2"
    paste "${OUT}/$(basename ${input}).words" "${OUT}/$(basename ${input}).nextwords" "${OUT}/$(basename ${input}).nextwords2" |
    sort | uniq -c  > "${OUT}/$(basename ${input}).trigrams"

done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
