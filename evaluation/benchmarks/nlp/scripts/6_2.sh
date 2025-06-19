#!/bin/bash
# tag: four-letter words
# set -e

# the original script has both versions
mkdir -p "$OUT"

for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -c '^....$' > ${OUT}${input}.stdout0
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | grep -c '^....$'  > ${OUT}${input}.stdout1
done

echo 'done';
# rm -rf "$OUT"
