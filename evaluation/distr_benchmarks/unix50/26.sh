#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN94=$IN_PRE/9.4.txt
# 9.4: four corners with E centered, for an "X" configuration
hdfs dfs -cat -ignoreCrc $IN94 | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'

