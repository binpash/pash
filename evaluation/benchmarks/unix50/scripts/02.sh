#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN1=$IN_PRE/1.txt
# 1.1: extract names and sort
cat "$IN1" | cut -d ' ' -f 2 | sort
