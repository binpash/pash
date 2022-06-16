#!/bin/bash
# tag: count_vowel_seq
# set -e 

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/2_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat -ignoreCrc $IN/$input | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c  > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
