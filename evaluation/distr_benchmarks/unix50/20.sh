#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN8=$IN_PRE/8.txt
# 8.3: find names of the four people most involved with unix
hdfs dfs -cat -ignoreCrc $IN8 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

