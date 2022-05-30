#!/bin/bash 
# tag: sort_words_by_num_of_syllables
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/8.1/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

pure_func() {
    input=$1
    cat > ${OUT}/${input}.words
    tr -sc '[AEIOUaeiou\012]' ' ' < ${OUT}/${input}.words | awk '{print NF}' > ${OUT}/${input}.syl
    paste ${OUT}/${input}.syl ${OUT}/${input}.words | sort -nr | sed 5q
}
export -f pure_func
for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | pure_func $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
