#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN11=$IN_PRE/11.txt
# 11.1: year Ritchie and Thompson receive the Hamming medal
hdfs dfs -cat $IN11 | grep 'UNIX' | cut -f 1

