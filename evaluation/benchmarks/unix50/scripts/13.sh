#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN5=$IN_PRE/05.txt
# 5.1: extract hello world
cat "$IN5" | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12
