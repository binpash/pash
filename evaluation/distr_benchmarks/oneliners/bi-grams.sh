#!/bin/bash
# Find all 2-grams in a piece of text

IN=${IN:-/1G.txt}

. bi-gram.aux.sh

hdfs dfs -cat $IN |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq


