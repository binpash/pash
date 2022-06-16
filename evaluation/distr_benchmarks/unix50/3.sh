#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN1=$IN_PRE/1.txt
# 1.2: extract names and sort
hdfs dfs -cat -ignoreCrc $IN1 | head -n 2 | cut -d ' ' -f 2

