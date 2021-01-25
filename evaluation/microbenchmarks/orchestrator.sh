#!/bin/bash

microbenchmark="comm-par-test2.sh diff.sh double_sort.sh for_loop_simple.sh fun-def.sh grep.sh micro_1000.sh minimal_grep.sh minimal_sort.sh no_in_script.sh sed-test.sh set-diff.sh shortest_scipts.sh sort.sh topn.sh wf.sh"
DIGITS=3
TIMEFORMAT="RadhasTimer %${DIGITS}R"

echo "Start" >> ./timings.txt

#            1     2     3   4   5  6  7 
# diff.sh  10.33  6.1  4.4 ....
# sort.sh  8.22   5.1 ...
# ...
echo '' > ./timings.csv

PER_SCRIPT="widths:"
for width in $(seq 8); do
  echo $width > perrun.csv
  PER_SCRIPT="$PER_SCRIPT,  $(cat perrun.csv)"
done
echo $PER_SCRIPT >> ./timings.csv

# 1 2 4 8
for f in $microbenchmark; do
  echo $f
  PER_SCRIPT="$f:"
  for width in $(seq 8); do
    echo $width
    (time $PASH_TOP/pa.sh -w $width $f) 2>&1 | grep RadhasTimer | sed 's/RadhasTimer //' > perrun.csv
    PER_SCRIPT="$PER_SCRIPT,  $(cat perrun.csv)"
  done
  echo $PER_SCRIPT >> ./timings.csv
done