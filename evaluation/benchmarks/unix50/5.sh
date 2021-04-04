#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN2=$IN_PRE/2.txt
# 2.1: get all Unix utilities
cat $IN2 | cut -d ' ' -f 4 | tr -d ','

