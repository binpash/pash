#!/bin/bash 
# tag: find_anagrams.sh
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
# need to generate words
cat ${IN}* | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
# need to generate types
sort -u ${OUT}.words > ${OUT}.types
# Actual find anagram script
rev < ${OUT}.types > ${OUT}.types.rev
sort ${OUT}.types ${OUT}.types.rev |
uniq -c |
awk "\$1 >= 2 {print \$2}"
