#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

IN=${IN:-/rep3_10M.txt}

hdfs dfs -cat -ignoreCrc $IN | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | tr A-Z a-z | sort | uniq -c | sort -rn 
