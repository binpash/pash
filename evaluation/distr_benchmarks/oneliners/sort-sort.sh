#!/bin/bash
# Calculate sort twice

IN=${IN:-/oneliners/1G.txt}

hdfs dfs -cat $IN | tr A-Z a-z | sort | sort -r
