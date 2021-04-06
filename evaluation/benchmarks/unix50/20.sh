#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN8=$IN_PRE/8.txt
# 8.3: find names of the four people most involved with unix
cat $IN8 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

