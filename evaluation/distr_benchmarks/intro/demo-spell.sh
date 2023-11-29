#!/bin/sh

cd "$(dirname $0)"

[ -z $PASH_TOP ] && { 
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}
DICT="$DISH_TOP/evaluation/distr_benchmarks/intro/input/sorted_words"
IN=${IN:-/intro/100M.txt}
hdfs dfs -cat -ignoreCrc $IN |
    tr A-Z a-z |
    tr -cs A-Za-z '\n' |
    sort |
    uniq |
    comm -13 $DICT -
