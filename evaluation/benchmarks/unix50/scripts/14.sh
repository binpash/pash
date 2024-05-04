#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/inputs}
IN6=$IN_PRE/06.txt
# 6.1: order the bodies by how easy it would be to land on them in Thompson's Space Travel game when playing at the highest simulation scale
cat "$IN6" | awk "{print \$2, \$0}" | sort -nr | cut -d ' ' -f 2
