#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.3: find pieces captured by Belle with a pawn
hdfs dfs -cat -ignoreCrc $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l

