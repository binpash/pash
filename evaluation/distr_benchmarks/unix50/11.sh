#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.5: 4.4 + pawns
hdfs dfs -cat -ignoreCrc $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr

