#!/bin/bash 
# tag: sort_words_by_num_of_syllables
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/8.1/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

run_tests() {
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}.words
    tr -sc '[AEIOUaeiou\012]' ' ' < ${OUT}/${input}.words | awk '{print NF}' > ${OUT}/${input}.syl
    paste ${OUT}/${input}.syl ${OUT}/${input}.words | sort -nr | sed 5q
}
export -f run_tests
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    run_tests $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
