#!/bin/bash
# Find all 2-grams in a piece of text

IN=${IN:-/oneliners/1G.txt}

. bi-gram.aux.sh

hdfs dfs -cat -ignoreCrc $IN |
  tr -c 'A-Za-z' '[\n*]' | 
  grep -v "^\s*$" |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq


