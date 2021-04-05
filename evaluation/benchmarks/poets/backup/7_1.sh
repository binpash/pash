#!/bin/bash
# tag: count_morphs
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
# FIXME: `spell -v' does not exist; instead of spell: sed 's/ly$/-ly/g'
ls $IN | xargs cat | spell | sed 's/ .*//g' | sort | uniq -c 
