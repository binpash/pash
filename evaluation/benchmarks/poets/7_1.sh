#!/bin/bash
# tag: count_morphs
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
# FIXME: `spell -v' does not exist; instead of spell: sed 's/ly$/-ly/g'
#ls ${IN} | sed "s;^;$IN;"| xargs cat | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c 
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/7_1}
mkdir -p ${OUT}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > ${OUT}/${input}.out
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm ${OUT}
