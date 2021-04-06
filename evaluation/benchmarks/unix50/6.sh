#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN3=$IN_PRE/3.txt
# 3.1: get lowercase first letter of last names (awk)
cat $IN3 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

