#!/bin/bash
# tag: sort_words_by_rhyming.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev 
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/3_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev > ${OUT}/${input}.out
done

#for output in $(ls ${OUT} | sed "s;^;$OUT;")
#do
#    cat $output
#done

echo 'done';
rm -rf "${OUT}"
