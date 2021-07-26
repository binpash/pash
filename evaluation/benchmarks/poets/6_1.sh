#!/bin/bash
# tag: trigram_rec
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}

trigrams() {
    tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
    tail +2 ${OUT}.words > ${OUT}.nextwords
    tail +3 ${OUT}.words > ${OUT}.nextwords2
    paste ${OUT}.words ${OUT}.nextwords ${OUT}.nextwords2 | sort | uniq -c
    rm -f ${OUT}.words ${OUT}.nextwords ${OUT}.nextwords2
}

ls ${IN} | sed "s;^;$IN;" | xargs cat | grep 'the land of' | trigrams | sort -nr | sed 5q
ls ${IN} | sed "s;^;$IN;" | xargs cat | grep 'And he said' | trigrams | sort -nr | sed 5q
