#!/bin/bash
# tag: four-letter words
# set -e

# the original script has both versions
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/6_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' 
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$' 


for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' > ${OUT}/${input}.out0
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$'  > ${OUT}/${input}.out1
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
