#!/bin/bash
# tag: compare_exodus_genesis.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/benchmarks/nlp/input/exodus}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/8.3_3/}
ENTRIES=${ENTRIES:-1060}
mkdir -p $OUT

run_tests() {
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}1.types
    tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > ${OUT}/${input}2.types
    sort $OUT/${input}1.types ${OUT}/${input}2.types ${OUT}/${input}2.types | uniq -c | head 

}
export -f run_tests
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    run_tests $input  > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
