#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN8=$IN_PRE/8.txt
# 8.4: find longest words without hyphens
hdfs dfs -cat $IN8 | tr -c "[a-z][A-Z]" '\n' | sort | awk "length >= 16"

