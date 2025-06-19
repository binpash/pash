#!/bin/bash
# tag: count_vowel_seq
# set -e

mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c  > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
