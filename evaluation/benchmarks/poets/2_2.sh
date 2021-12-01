#!/bin/bash
# tag: count_vowel_seq
# set -e 

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c 
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/2_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
echo $ENTRIES
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c  > ${OUT}/${input}.out
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
