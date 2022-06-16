#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.6: piece used the most by Belle
hdfs dfs -cat -ignoreCrc $IN4 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

