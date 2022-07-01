#!/bin/bash
export IN_PRE=${IN_PRE:-/unix50}
IN4=$IN_PRE/4.txt
# 4.4: histogram of Belle's captures (-pawns) by each type of piece
hdfs dfs -cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr

