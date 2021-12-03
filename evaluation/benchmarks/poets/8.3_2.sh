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

OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.3_2/}
mkdir -p "$OUT"
ENTRIES=${ENTRIES:-1000}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.words
    sort -u ${OUT}/${input}.words > ${OUT}/${input}.types
    rev < ${OUT}/${input}.types > ${OUT}/${input}.types.rev
    sort ${OUT}/${input}.types ${OUT}/${input}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}" > ${OUT}/${input}.out
done

#for output in $(ls ${OUT} | sed "s;^;$OUT;")
#do
#    cat $output
#done

echo 'done';
rm -rf "${OUT}"
