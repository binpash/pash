#!/bin/bash 
# tag: sort_words_by_num_of_syllables
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}

ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}.words
tr -sc '[AEIOUaeiou\012]' ' ' < ${OUT}.words | awk '{print NF}' > ${OUT}.syl
paste ${OUT}.syl ${OUT}.words | sort -nr | sed 5q
