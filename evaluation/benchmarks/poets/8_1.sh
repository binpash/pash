#!/bin/bash 
# tag: sort_words_by_num_of_syllables
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}.words
#tr -sc '[AEIOUaeiou\012]' ' ' < ${OUT}.words | awk '{print NF}' > ${OUT}.syl
#paste ${OUT}.syl ${OUT}.words | sort -nr | sed 5q
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.1/}
mkdir -p "$OUT"
ENTRIES=${ENTRIES:-1000}
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}.words
    tr -sc '[AEIOUaeiou\012]' ' ' < ${OUT}/${input}.words | awk '{print NF}' > ${OUT}/${input}.syl
    paste ${OUT}/${input}.syl ${OUT}/${input}.words | sort -nr | sed 5q > ${OUT}/${input}.out
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
