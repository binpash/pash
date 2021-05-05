#!/bin/bash 
# tag: bigrams.sh
# set -e

# Bigrams (contrary to our version, this uses intermediary files)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.input.words
tail +2 ${OUT}.input.words > ${OUT}.input.nextwords
paste ${OUT}.input.words ${OUT}.input.nextwords | sort | uniq -c > ${OUT}.input.bigrams
