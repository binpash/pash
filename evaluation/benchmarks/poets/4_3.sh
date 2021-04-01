#!/bin/bash 
# tag: bigrams.sh
set -e

# Bigrams (contrary to our version, this uses intermediary files)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > input.words
tail +2 input.words > input.nextwords
paste input.words input.nextwords | sort | uniq -c > input.bigrams
