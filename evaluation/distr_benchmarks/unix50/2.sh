#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN1=$IN_PRE/1.txt
# 1.1: extract names and sort
hdfs dfs -cat $IN1 | cut -d ' ' -f 2 | sort

