#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN96=$IN_PRE/9.6.txt
# 9.6: Follow the directions for grep
cat $IN96 | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'

