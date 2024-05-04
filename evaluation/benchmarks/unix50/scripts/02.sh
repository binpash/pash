#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN1=$IN_PRE/01.txt
# 1.1: extract names and sort
cat "$IN1" | cut -d ' ' -f 2 | sort
