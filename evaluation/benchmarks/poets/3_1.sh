#!/bin/bash
# tag: sort
# set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat| tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | sort -nr
