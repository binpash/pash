#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN9_6=$IN_PRE/09.6.txt
# 9.6: Follow the directions for grep
cat "$IN9_6" | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'
