#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN3=$IN_PRE/3.txt
# 3.1: get lowercase first letter of last names (awk)
hdfs dfs -cat -ignoreCrc $IN3 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

