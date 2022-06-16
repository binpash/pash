#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.2: find pieces captured by Belle
hdfs dfs -cat -ignoreCrc $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

