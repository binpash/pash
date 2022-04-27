#!/bin/bash
# Sort input

IN=${IN:-/1G.txt}

hdfs dfs -cat $IN | sort

