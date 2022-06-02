#!/bin/bash 
# tag: bigrams_appear_twice.sh
# set -e

# Calculate the bigrams (based on 4_3.sh script)
IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/8.2_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.input.words
    tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
    paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords | sort | uniq -c > ${TEMPDIR}/${input}.input.bigrams
    awk "\$1 == 2 {print \$2, \$3}" ${TEMPDIR}/${input}.input.bigrams
    rm -rf {TEMPDIR}
}

export -f pure_func
for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
