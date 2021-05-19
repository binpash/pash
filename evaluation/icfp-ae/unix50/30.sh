#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN98=$IN_PRE/9.8.txt
# 9.8: TELE-communications
cat $IN98 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'

