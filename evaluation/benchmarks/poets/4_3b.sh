#!/bin/bash 
#tag: count_trigrams.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
tail +2 ${OUT}.words > ${OUT}.nextwords
tail +3 ${OUT}.words > ${OUT}.nextwords2
paste ${OUT}.words ${OUT}.nextwords ${OUT}.nextwords2 |
sort | uniq -c  > ${OUT}.trigrams
