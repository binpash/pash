#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN8=$IN_PRE/8.txt
# 8.1: count unix birth-year
cat $IN8 | tr ' ' '\n' | grep 1969 | wc -l

