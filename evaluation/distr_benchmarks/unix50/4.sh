#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN1=$IN_PRE/1.txt
# 1.3: sort top first names
hdfs dfs -cat $IN1 | cut -d ' ' -f 1 | sort | uniq -c | sort -r

