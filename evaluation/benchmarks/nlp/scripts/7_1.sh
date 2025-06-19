#!/bin/bash
# tag: count_morphs
# set -e
mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf ${OUT}
