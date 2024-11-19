#!/bin/bash
# tag: trigram_rec
# set -e


pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    tr -sc '[A-Z][a-z]' '[\012*]' > ${TEMPDIR}/${input}.words
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords
    tail +3 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords2
    paste ${TEMPDIR}/${input}.words ${TEMPDIR}/${input}.nextwords ${TEMPDIR}/${input}.nextwords2 | sort | uniq -c
    rm -rf ${TEMPDIR}
}
export -f pure_func

for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | grep 'the land of' | pure_func ${input} | sort -nr | sed 5q > ${OUT}${input}.0.out
    cat $IN$input | grep 'And he said' | pure_func ${input} | sort -nr | sed 5q > ${OUT}${input}.1.out
done

echo 'done';
# rm -rf "$OUT"
