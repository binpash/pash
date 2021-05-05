#!/bin/bash
# tag: merge_upper
# set -e

# Merge upper and lower counts
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c
