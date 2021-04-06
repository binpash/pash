#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN8=$IN_PRE/8.txt
# 8.2: find Bell Labs location where Dennis Ritchie had his office
cat $IN8 | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk "{\$1=\$1};1"

