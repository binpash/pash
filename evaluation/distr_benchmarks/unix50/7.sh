#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.1: find number of rounds
hdfs dfs -cat $IN4 | tr ' ' '\n' | grep '\.' | wc -l

