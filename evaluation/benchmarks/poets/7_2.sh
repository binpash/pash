#!/bin/bash
# set -e
# tag: count_consonant_sequences

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr '[a-z]' '[A-Z]' | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c
