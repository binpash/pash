#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN7=$IN_PRE/7.txt
# 7.1: identify number of AT&T unix versions
hdfs dfs -cat -ignoreCrc $IN7 | cut -f 1 | grep 'AT&T' | wc -l

