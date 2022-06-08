#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN99=$IN_PRE/9.9.txt
# 9.9:
hdfs dfs -cat $IN99 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'

