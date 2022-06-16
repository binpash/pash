#!/bin/bash
# tag: 2-syllable words
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/6_5/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat -ignoreCrc $IN/$input  | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -i '^[^aeiou]*[aeiou][^aeiou]*[aeiou][^aeiou]$' | sort | uniq -c | sed 5q > ${OUT}${input}.out
done

echo 'done';
rm -rf "${OUT}"
