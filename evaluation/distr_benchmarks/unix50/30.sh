#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN98=$IN_PRE/9.8.txt
# 9.8: TELE-communications
hdfs dfs -cat -ignoreCrc $IN98 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'

