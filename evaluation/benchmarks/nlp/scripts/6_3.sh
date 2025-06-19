#!/bin/bash
# tag: words_no_vowels
# set -e

mkdir -p "$OUT"

for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -vi '[aeiou]' | sort | uniq -c > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
