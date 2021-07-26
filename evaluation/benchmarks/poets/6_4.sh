#!/bin/bash
# tag: 1-syllable words
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[ 12*]' | grep -i '^[^aeiou]*[aeiou][^aeiou]*$' | sort | uniq -c | sed 5q
