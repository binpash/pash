#!/bin/bash

# tag: 2-syllable words
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
cat ${IN}* | tr -sc '[A-Z][a-z]' ' [\012*]' | grep -i '^[^aeiou]*[aeiou][^aeiou]*[aeiou][^aeiou]$' | sort | uniq -c | sed 5q
