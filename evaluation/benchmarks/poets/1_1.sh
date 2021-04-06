#!/bin/bash
# tag: count_words
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
cat ${IN}* | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c
