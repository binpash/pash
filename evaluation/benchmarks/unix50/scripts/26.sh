#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN9_4=$IN_PRE/9.4.txt
# 9.4: four corners with E centered, for an "X" configuration
cat "$IN9_4" | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'
