#!/bin/bash
# tag: compare_exodus_genesis.sh
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/benchmarks/poets/input/exodus}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u >  1.types
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > 2.types
# FIXME do we really need the same thing twice?
sort 1.types 2.types 2.types | uniq -c | head
