#!/bin/bash
# Sort input

IN=${IN:-/oneliners/1G.txt}

hdfs dfs -cat $IN | sort

