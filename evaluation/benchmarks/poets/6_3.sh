#!/bin/bash
# tag: words_no_vowels
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c
