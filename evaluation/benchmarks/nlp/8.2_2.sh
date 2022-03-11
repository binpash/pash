#!/bin/bash 
# tag: bigrams_appear_twice.sh
# set -e

# Calculate the bigrams (based on 4_3.sh script)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/8.2_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

run_tests() {
    input=$1
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.input.words
    tail +2 ${OUT}/${input}.input.words > ${OUT}/${input}.input.nextwords
    paste ${OUT}/${input}.input.words ${OUT}/${input}.input.nextwords | sort | uniq -c > ${OUT}/${input}.input.bigrams
    awk "\$1 == 2 {print \$2, \$3}" ${OUT}/${input}.input.bigrams
}

export -f run_tests
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    run_tests $input  > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
