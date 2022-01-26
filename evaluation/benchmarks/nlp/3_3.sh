#!/bin/bash
# tag: sort_words_by_rhyming.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/3_3/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev > ${OUT}/${input}.out
done

echo 'done';
rm -rf "${OUT}"
