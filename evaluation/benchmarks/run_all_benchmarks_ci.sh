#!/bin/bash

## This script is necessary to ensure that sourcing happens with bash
source run.seq.sh
source run.par.sh

compare_outputs(){
  dir=$1
  outputs=$(ls $dir | grep "seq" | sed 's/.seq.out$//')
  for out in $outputs;
  do
    seq_output="${dir}/${out}.seq.out"
    pash_output="${dir}/${out}.par.out"
    diff -q "$seq_output" "$pash_output"
  done
}

#oneliners
#oneliners_pash
#
#compare_outputs "oneliners/outputs"
#
#unix50
#unix50_pash
#
#compare_outputs "unix50/outputs"
#
#poets
#poets_pash
#
#compare_outputs "poets/outputs"
#
#web-index
#web-index_pash
#
#compare_outputs "web-index/outputs"
rm -f $1/*.res
echo "group,Bash,Pash16" > results.time
b=$($1      | grep -v executing | sed 's\.sh:\_bash\g' |  sed 's\,\.\g' | awk '{ print $1 "," $2}')
p=$($1_pash | grep -v executing | sed 's\.sh:\_pash\g' |  sed 's\,\.\g' | awk '{ print $1 "," $2}' | awk '{sub(/[^,]*/,"");sub(/,/,"")} 1')
res=$(paste -d '@' <(echo "$b") <(echo "$p"))
echo "$res" | sed 's\@\,\g' >> results.time
compare_outputs "$1/outputs"
cat results.time
