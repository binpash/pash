#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN9_8=$IN_PRE/09.8.txt
# 9.8: TELE-communications
cat "$IN9_8" | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'
