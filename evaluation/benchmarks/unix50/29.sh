#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN97=$IN_PRE/9.7.txt
# 9.7: Four corners
cat $IN97 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

