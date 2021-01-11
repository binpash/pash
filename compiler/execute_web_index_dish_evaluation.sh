#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}


eval_dir="../evaluation/"
directory="${eval_dir}/scripts/web-index/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"

output_dir="/tmp/web-index-output"
p1_out="${output_dir}/p1.out"

## TODO: Move the scripts outside of the intermediary
intermediary_dir="${eval_dir}/intermediary/"
results_dir="${eval_dir}/results/"

input_dir="${HOME}/wikipedia/"

## Make the temporary output dir
mkdir -p $output_dir

# cp "$input_dir/index_h_100000.txt" $p1_out
cp "$input_dir/index_h_100.txt" $p1_out


## TODO: If we have time, we need to remove this split

## Split the input file
split -n l/2 -d $p1_out ${p1_out}_2_
split -n l/16 -d $p1_out ${p1_out}_16_

echo "Intermediate files have been successfully produced."

export IN_DIR="${output_dir}"
export WIKI="${input_dir}"
export WEB_INDEX_DIR="${directory}"

web_index_2_script="${intermediary_dir}/web-index_full_2_seq.sh"
web_index_16_script="${intermediary_dir}/web-index_full_16_seq.sh"

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
{ time /bin/bash $web_index_2_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

echo "Executing the script with pash -w 2 (log in: /tmp/pash_2_log)"
{ time $PASH_TOP/pa.sh -w 2 --log_file /tmp/pash_2_log --output_time $web_index_2_script ; } 1> "$pash_width_2_output" 2> >(tee "${pash_width_2_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_2_output | head

echo "Executing the script with pash -w 16 (log in: /tmp/pash_16_log)"
{ time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $web_index_16_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
