#!/bin/bash

dir="../scripts/max-temp/"
results="../evaluation/results/"
# temp_dir="/dev/shm/max-temp/"
temp_dir="./temp/"

mkdir -p $temp_dir
time {
seq 2005 2005 |
  sed 's;^;http://ndr.md/noaa/;' |
  sed 's;$;/;' |
  xargs -n 1 curl -s > "$temp_dir/p1.txt" ; } 2> >(tee "${results}max-temp-p1.time" >&2)

time {
cat "$temp_dir/p1.txt" |
  grep gz |
  tr -s ' ' |
  cut -d ' ' -f9 > "$temp_dir/p2.txt" ; } 2> >(tee "${results}max-temp-p2.time" >&2)

time {
cat "$temp_dir/p2.txt" |
    sed 's;^;http://ndr.md/noaa/2005/;' > "$temp_dir/p3a.txt" ; } 2> >(tee "${results}max-temp-p3a.time" >&2)

distr_output_dir=/tmp/distr_output

./execute_compile_evaluation_script.sh "p3b_1" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "p3b_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "p3b_10" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "p3b4_1" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "p3b4_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "p3b4_10" $distr_output_dir 1 0

