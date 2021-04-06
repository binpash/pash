#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN4=$IN_PRE/4.txt
# 4.2: find pieces captured by Belle
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

