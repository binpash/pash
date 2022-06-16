#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN93=$IN_PRE/9.3.txt
# 9.3: animal that used to decorate the Unix room
hdfs dfs -cat -ignoreCrc $IN93 | cut -c 1-2 | tr -d '\n'

