#!/bin/bash
# tag: merge_upper
# set -e

# Merge upper and lower counts
mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr '[a-z]' '[A-Z]' |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
