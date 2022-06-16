#!/bin/bash
# Match complex regular-expression over input

IN=${IN:-/oneliners/1G.txt}

hdfs dfs -cat -ignoreCrc $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
