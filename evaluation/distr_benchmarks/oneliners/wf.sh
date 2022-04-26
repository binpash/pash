#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

IN=${IN:-/10M.txt}

hdfs dfs -cat $IN |  tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn 
