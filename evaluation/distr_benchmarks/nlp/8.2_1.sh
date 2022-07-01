#!/bin/bash
# tag: vowel_sequences_gr_1K.sh
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/8.2_1/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000" > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
