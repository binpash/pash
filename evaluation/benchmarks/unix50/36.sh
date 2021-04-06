#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN11=$IN_PRE/11.txt
# 11.2: most repeated first name in the list?
cat $IN11 | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d

