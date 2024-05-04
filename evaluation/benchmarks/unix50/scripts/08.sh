#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN4=$IN_PRE/04.txt
# 4.2: find pieces captured by Belle
cat "$IN4" | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l
