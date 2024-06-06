#!/bin/bash
# tag: sort_words_by_rhyming.sh
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/3_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c | rev | sort | rev > ${OUT}/${input}.out
done

echo 'done';
# rm -rf "$OUT"
