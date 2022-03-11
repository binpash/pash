#!/bin/bash
# tag: trigram_rec
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/6_1/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

trigrams() {
    input=$1
    tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.words
    tail +2 ${OUT}/${input}.words > ${OUT}/${input}.nextwords
    tail +3 ${OUT}/${input}.words > ${OUT}/${input}.nextwords2
    paste ${OUT}/${input}.words ${OUT}/${input}.nextwords ${OUT}/${input}.nextwords2 | sort | uniq -c
    rm -f ${OUT}/${input}.words ${OUT}/${input}.nextwords ${OUT}/${input}.nextwords2
}
export -f trigrams

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN"/"$input | grep 'the land of' | trigrams $IN/${input} | sort -nr | sed 5q > ${OUT}/${input}.out0
    cat $IN"/"$input | grep 'And he said' | trigrams $IN/${input} | sort -nr | sed 5q > ${OUT}/${input}.out1
done

echo 'done';
rm -rf "${OUT}"
