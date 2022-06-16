#!/bin/bash 
#tag: count_trigrams.sh
# set -e

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/4_3b/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.words
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords2
    paste ${TEMPDIR}/${input}.words ${TEMPDIR}/${input}.nextwords ${TEMPDIR}/${input}.nextwords2 |
    sort | uniq -c 
    rm -rf ${TEMPDIR}
}
export -f pure_func
for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat -ignoreCrc $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}/${input}.trigrams
done

echo 'done';
rm -rf ${OUT}
