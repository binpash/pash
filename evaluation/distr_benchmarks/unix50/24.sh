#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN92=$IN_PRE/9.2.txt
# 9.2: extract the word BELL
hdfs dfs -cat $IN92 | cut -c 1-1 | tr -d '\n'

