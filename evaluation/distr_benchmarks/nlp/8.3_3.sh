#!/bin/bash
# tag: compare_exodus_genesis.sh
# set -e

IN=${IN:-/nlp/pg/}
INPUT2=${INPUT2:-/nlp/exodus}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/8.3_3/}
ENTRIES=${ENTRIES:-1060}
mkdir -p $OUT

pure_func() {
    input=$1
    cat > ${OUT}/${input}1.types
    hdfs dfs -cat  ${INPUT2} | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}2.types
    sort $OUT/${input}1.types ${OUT}/${input}2.types ${OUT}/${input}2.types | uniq -c | head 

}
export -f pure_func
for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | pure_func $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
