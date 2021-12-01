#!/bin/bash 
# tag: bigrams.sh
# set -e

# Bigrams (contrary to our version, this uses intermediary files)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/4_3}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.input.words
#tail +2 ${OUT}.input.words > ${OUT}.input.nextwords
#paste ${OUT}.input.words ${OUT}.input.nextwords | sort | uniq -c > ${OUT}.input.bigrams
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.input.words
    tail +2  ${OUT}/${input}.input.words > ${OUT}/${input}.input.nextwords
    paste ${OUT}/${input}.input.words ${OUT}/${input}.input.nextwords | sort | uniq -c > ${OUT}/${input}.input.bigrams
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm -rf ${OUT}
