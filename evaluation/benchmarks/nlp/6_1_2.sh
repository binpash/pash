#!/bin/bash
# tag: uppercase_by_type
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/nlp/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/nlp/output/6_1_2/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^[A-Z]' > ${OUT}/${input}.out
done

echo 'done';
rm -rf ${OUT}
