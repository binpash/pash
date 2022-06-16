#!/bin/bash
# tag: count_morphs
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/7_1/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat -ignoreCrc $IN/$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
rm -rf ${OUT}
