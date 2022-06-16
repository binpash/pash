#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN8=$IN_PRE/8.txt
# 8.1: count unix birth-year
hdfs dfs -cat -ignoreCrc $IN8 | tr ' ' '\n' | grep 1969 | wc -l

