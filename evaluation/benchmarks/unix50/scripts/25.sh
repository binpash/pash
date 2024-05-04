#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN9_3=$IN_PRE/09.3.txt
# 9.3: animal that used to decorate the Unix room
cat "$IN9_3" | cut -c 1-2 | tr -d '\n'
