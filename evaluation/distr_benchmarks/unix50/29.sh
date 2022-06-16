#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN97=$IN_PRE/9.7.txt
# 9.7: Four corners
hdfs dfs -cat -ignoreCrc $IN97 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

