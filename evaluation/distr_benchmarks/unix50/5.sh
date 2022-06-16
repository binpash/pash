#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN2=$IN_PRE/2.txt
# 2.1: get all Unix utilities
hdfs dfs -cat -ignoreCrc $IN2 | cut -d ' ' -f 4 | tr -d ','

