#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN91=$IN_PRE/9.1.txt
# 9.1: extract the word PORT
hdfs dfs -cat $IN91 | tr ' ' '\n' | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4

