#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN10=$IN_PRE/10.txt
# 10.3: extract Ritchie's username
cat $IN10 | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

