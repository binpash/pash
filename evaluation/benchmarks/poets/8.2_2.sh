#!/bin/bash 
# tag: bigrams_appear_twice.sh
# set -e

# Calculate the bigrams (based on 4_3.sh script)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.input.words
tail +2 ${OUT}.input.words > ${OUT}.input.nextwords
paste ${OUT}.input.words ${OUT}.input.nextwords | sort | uniq -c > ${OUT}.input.bigrams
# find the bigrams that appear exactly twice 
awk "\$1 == 2 {print \$2, \$3}" ${OUT}.input.bigrams
