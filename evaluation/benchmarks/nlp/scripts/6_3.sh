#!/bin/bash
# tag: words_no_vowels
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -vi '[aeiou]' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
# rm -rf "$OUT"
