#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN6=$IN_PRE/6.txt
# 6.1: order the bodies by how easy it would be to land on them in Thompson's Space Travel game when playing at the highest simulation scale
hdfs dfs -cat -ignoreCrc $IN6 | awk "{print \$2, \$0}" | sort -nr | cut -d ' ' -f 2

