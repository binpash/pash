#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN8=$IN_PRE/8.txt
# 8.2: find Bell Labs location where Dennis Ritchie had his office
hdfs dfs -cat -ignoreCrc $IN8 | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk "{\$1=\$1};1"

