#!/bin/bash
# tag: 1-syllable words
# set -e

mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat ${IN}/${input} | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -i '^[^aeiou]*[aeiou][^aeiou]*$' | sort | uniq -c | sed 5q > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
