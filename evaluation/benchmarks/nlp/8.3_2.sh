#!/bin/bash 
# tag: find_anagrams.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/8.3_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

run_tests() {
    input=$1
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${OUT}/${input}.types
    rev < ${OUT}/${input}.types > ${OUT}/${input}.types.rev
    sort ${OUT}/${input}.types ${OUT}/${input}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}"
}

export -f run_tests
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    run_tests $input > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
