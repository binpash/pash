#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN9_2=$IN_PRE/9.2.txt
# 9.2: extract the word BELL
cat "$IN9_2" | cut -c 1-1 | tr -d '\n'
