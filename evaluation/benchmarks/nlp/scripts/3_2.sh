#!/bin/bash
# tag: sort_words_by_folding
# set -e

mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c | sort -f > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf ${OUT}
