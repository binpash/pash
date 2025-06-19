#!/bin/bash
# tag: uppercase_by_type
# set -e
mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | grep -c '^[A-Z]' > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf ${OUT}
