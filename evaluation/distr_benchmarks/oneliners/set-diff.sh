#!/bin/bash
# Show the set-difference between two streams (i.e., elements in the first that are not in the second).
# https://stackoverflow.com/questions/2509533/bash-linux-set-difference-between-two-text-files

IN=${IN:-/oneliners/1G.txt}

mkfifo s1 s2

hdfs dfs -cat $IN |
    cut -d ' ' -f 1 |
    tr [:lower:] [:upper:] |
    sort > s1 &

hdfs dfs -cat $IN |
    cut -d ' ' -f 1 |
    sort > s2 &

comm -23 s1 s2

rm s1 s2
