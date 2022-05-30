#!/bin/bash
# tag: sort_words_by_folding
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/3_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c | sort -f > ${OUT}/${input}
done

echo 'done';
rm -rf ${OUT}
