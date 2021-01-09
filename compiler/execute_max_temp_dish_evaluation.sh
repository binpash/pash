#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

eval_dir="../evaluation/"
results_dir="${eval_dir}/results/"

start_year=2000
# end_year=2004 # For the full evaluation
end_year=2000

max_temp_complete_script="${eval_dir}/scripts/max-temp-complete.sh"

seq_output=/tmp/seq_output
pash_width_16_output=/tmp/pash_16_output
seq_time="$results_dir/max-temp-complete-seq.time"
pash_width_16_time="$results_dir/max-temp-complete-16-pash.time"

echo "Executing the complete max-temp script with bash..."
seq "$start_year" "$end_year" | { time /bin/bash $max_temp_complete_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

## TODO: Figure a way to redirect stderr from curl not having input (maybe use an argument for time)
echo "Executing the complete max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_complete_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head

## TODO: Add an optional flag for the following
echo "Executing preprocessing and processing separately"

max_temp_preprocess_script="${eval_dir}/scripts/max-temp-preprocess.sh"
max_temp_process_script="${eval_dir}/scripts/max-temp-process.sh"

seq_preprocess_time="$results_dir/max-temp-preprocess-seq.time"
pash_width_16_preprocess_time="$results_dir/max-temp-preprocess-16-pash.time"
seq_process_time="$results_dir/max-temp-process-seq.time"
pash_width_16_process_time="$results_dir/max-temp-process-16-pash.time"
preprocess_output=/tmp/max-temp-preprocess-output

echo "Executing the preprocessing max-temp script with bash..."
seq "$start_year" "$end_year" | { time /bin/bash $max_temp_preprocess_script > $seq_output ; } 2> >(tee "${seq_preprocess_time}" >&2)

## TODO: Figure a way to redirect stderr from curl not having input
echo "Executing the preprocessing max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_preprocess_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_preprocess_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head

## Copy the sequential preprocess output to another file so that it doesn't get overwritten
cp "$seq_output" "${preprocess_output}"

## Export the input variable for the process script
export IN="${preprocess_output}"

echo "Executing the processing max-temp script with bash..."
{ time /bin/bash $max_temp_process_script > $seq_output ; } 2> >(tee "${seq_process_time}" >&2)

echo "Executing the processing max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
{ time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_process_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_process_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
