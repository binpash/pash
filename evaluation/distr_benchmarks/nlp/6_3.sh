#!/bin/bash
# tag: words_no_vowels
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/6_3/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
