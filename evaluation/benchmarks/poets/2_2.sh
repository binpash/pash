#!/bin/bash
# tag: count_vowel_seq
set -e 

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c 
