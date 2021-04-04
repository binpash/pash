#!/bin/bash
export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
IN7=$IN_PRE/7.txt
# 7.3: all the decades in which a unix version was released
cat $IN7 | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/

