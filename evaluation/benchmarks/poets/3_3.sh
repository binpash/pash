#!/bin/bash
# tag: sort_words_by_rhyming.sh
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
cat ${IN}* | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev 
