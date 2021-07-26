#!/bin/bash
# tag: uppercase_by_type
# set -e

# Uppercase words by type
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^[A-Z]'
