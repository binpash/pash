#!/bin/bash
# tag: count_words

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/1_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c > $1/${input}.out
done

echo 'done';
# rm -rf "$OUT"
