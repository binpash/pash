#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN11=$IN_PRE/11.txt
# 11.1: year Ritchie and Thompson receive the Hamming medal
cat $IN11 | grep 'UNIX' | cut -f 1

