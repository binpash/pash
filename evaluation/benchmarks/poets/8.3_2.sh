#!/bin/bash 
# tag: find_anagrams.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
## need to generate words
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
## need to generate types
#sort -u ${OUT}.words > ${OUT}.types
## Actual find anagram script
#rev < ${OUT}.types > ${OUT}.types.rev
#sort ${OUT}.types ${OUT}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}"

OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.3_2}
mkdir -p "$OUT"
ENTRIES=${ENTRIES:-1000}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000"
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr -sc '[A-Z][a-z]' '[\012*]' > "${OUT}/$(basename ${input}).words"
    sort -u "${OUT}/$(basename ${input}).words" > "${OUT}/$(basename ${input}).types"
    rev < "${OUT}/$(basename ${input}).types" > "${OUT}/$(basename ${input}).types.rev"
    sort "${OUT}/$(basename ${input}).types" "${OUT}/$(basename ${input}).types.rev" | uniq -c | awk "\$1 >= 2 {print \$2}" > "${OUT}/$(basename ${input}).out"
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
