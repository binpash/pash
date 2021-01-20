#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

eval_dir="../evaluation/"
directory="${eval_dir}/scripts/web-index/"
results_dir="${eval_dir}/results/"
input_dir="${HOME}/wikipedia/"

# export IN="$input_dir/index_h_100000.txt" # A couple hours in bash 
# export IN="$input_dir/index_h_1000.txt" # About 1.5 minutes in bash 
export IN="$input_dir/index_h_100.txt" # About 10 seconds in bash 
export WIKI="${input_dir}"
export WEB_INDEX_DIR="${directory}"

web_index_script="${eval_dir}/scripts/web-index.sh"

seq_output=/tmp/seq_output
pash_width_2_output=/tmp/pash_2_output
pash_width_16_output=/tmp/pash_16_output
seq_time="$results_dir/web-index-seq.time"
pash_width_2_time="$results_dir/web-index-2-pash.time"
pash_width_16_time="$results_dir/web-index-16-pash.time"

## Note: The sequential p2 stage takes about 7 seconds for 100! urls
##       This means that for the whole wikipedia (1mil urls) it will take about
##       20+ hours to complete.
echo "Executing the script with bash..."
{ time /bin/bash $web_index_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

echo "Executing the script with pash -w 2 (log in: /tmp/pash_2_log)"
{ time $PASH_TOP/pa.sh -w 2 --log_file /tmp/pash_2_log --output_time $web_index_script ; } 1> "$pash_width_2_output" 2> >(tee "${pash_width_2_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_2_output | head

echo "Executing the script with pash -w 16 (log in: /tmp/pash_16_log)"
{ time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $web_index_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
