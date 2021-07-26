#!/bin/bash
# tag: vowel_sequences_gr_1K.sh
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000"
