#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN4=$IN_PRE/4.txt
# 4.6: piece used the most by Belle
cat $IN4 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

