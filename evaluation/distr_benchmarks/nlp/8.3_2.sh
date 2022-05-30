#!/bin/bash 
# tag: find_anagrams.sh
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/8.3_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

pure_func() {
    input=$1
    sort -u > ${OUT}/${input}.types
    rev < ${OUT}/${input}.types > ${OUT}/${input}.types.rev
    sort ${OUT}/${input}.types ${OUT}/${input}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}"
}

export -f pure_func
for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
