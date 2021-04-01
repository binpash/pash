#!/bin/bash 
#tag: count_trigrams.sh
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${INPUT}.words
tail +2 ${INPUT}.words > ${INPUT}.nextwords
tail +3 ${INPUT}.words > ${INPUT}.nextwords2
paste ${INPUT}.words ${INPUT}.nextwords ${INPUT}.nextwords2 |
sort | uniq -c  > ${INPUT}.trigrams
