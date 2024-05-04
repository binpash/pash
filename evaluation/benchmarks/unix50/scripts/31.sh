#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN9_9=$IN_PRE/09.9.txt
# 9.9:
cat "$IN9_9" | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'
