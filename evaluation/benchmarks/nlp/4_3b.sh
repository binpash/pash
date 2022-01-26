#!/bin/bash 
#tag: count_trigrams.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/4_3b/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

run_tests() {
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.words
    tail +2 ${OUT}/${input}.words > ${OUT}/${input}.nextwords
    tail +2 ${OUT}/${input}.words > ${OUT}/${input}.nextwords2
    paste ${OUT}/${input}.words ${OUT}/${input}.nextwords ${OUT}/${input}.nextwords2 |
    sort | uniq -c 
}
export -f run_tests
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    run_tests $input > ${OUT}/${input}.trigrams
done

echo 'done';
rm -rf ${OUT}
