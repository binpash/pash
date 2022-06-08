#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN10=$IN_PRE/10.txt
# 10.1: count Turing award recipients while working at Bell Labs
hdfs dfs -cat $IN10 | sed 1d | grep 'Bell' | cut -f 2 | wc -l

