#!/bin/bash
# tag: count_morphs
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/7_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
# rm -rf ${OUT}
