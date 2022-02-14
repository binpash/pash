#!/bin/bash
# tag: four-letter words
# set -e

# the original script has both versions
IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/6_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' > ${OUT}/${input}.out0
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$'  > ${OUT}/${input}.out1
done

echo 'done';
rm -rf "${OUT}"
