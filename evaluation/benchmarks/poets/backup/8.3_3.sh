#!/bin/bash
# tag: compare_exodus_genesis.sh
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/benchmarks/poets/input/exodus}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u >  ${OUT}1.types
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > ${OUT}2.types
# FIXME do we really need the same thing twice?
sort ${OUT}1.types ${OUT}2.types ${OUT}2.types | uniq -c | head
