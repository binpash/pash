#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN1=$IN_PRE/01.txt
# 1.3: sort top first names
cat "$IN1" | cut -d ' ' -f 1 | sort | uniq -c | sort -r
