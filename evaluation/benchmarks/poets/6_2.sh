#!/bin/bash
# tag: four-letter words
set -e

# the original script has both versions
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' 
ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$' 

#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^....$' | paste -sd+ | bc
#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^....$' | paste -sd+ | bc

