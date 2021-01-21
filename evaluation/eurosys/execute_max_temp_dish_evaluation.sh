#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## There are two possible execution levels:
## options -s: end_year=2000
## options -l: end_year=2004 (The EuroSys evaluation)
start_year=2000
end_year=2000 # For the small evaluation

while getopts 'slh' opt; do
    case $opt in
        s) end_year=2000 ;;
        l) end_year=2004 ;;
        h) echo "There are three possible execution levels:"
           echo "option -s: end_year=2000"
           echo "option -l: end_year=2004 (The EuroSys evaluation)"
           exit 0 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"


eval_dir="$PASH_TOP/evaluation/"
results_dir="${eval_dir}/results/"


echo "Running max-temp evaluation for years: $start_year-$end_year"

max_temp_complete_script="${eval_dir}/scripts/max-temp-complete.sh"

seq_output=/tmp/seq_output
pash_width_16_output=/tmp/pash_16_output
seq_time="$results_dir/max-temp-complete-$start_year-$end_year-seq.time"
pash_width_16_time="$results_dir/max-temp-complete-$start_year-$end_year-16-pash.time"

echo "Executing the complete max-temp script with bash..."
seq "$start_year" "$end_year" | { time /bin/bash $max_temp_complete_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

## TODO: Figure a way to redirect stderr from curl not having input (maybe use an argument for time)
echo "Executing the complete max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_complete_script ; } 1> "$pash_width_16_output" 2> >(tee "$pash_width_16_time" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
 
## TODO: Add an optional flag for the following
echo "Executing preprocessing and processing separately"

max_temp_preprocess_script="${eval_dir}/scripts/max-temp-preprocess.sh"
max_temp_process_script="${eval_dir}/scripts/max-temp-process.sh"

seq_preprocess_time="$results_dir/max-temp-preprocess-$start_year-$end_year-seq.time"
pash_width_16_preprocess_time="$results_dir/max-temp-preprocess-$start_year-$end_year-16-pash.time"
seq_process_time="$results_dir/max-temp-process-$start_year-$end_year-seq.time"
pash_width_16_process_time="$results_dir/max-temp-process-$start_year-$end_year-16-pash.time"
preprocess_output=/tmp/max-temp-preprocess-output

echo "Executing the preprocessing max-temp script with bash..."
seq "$start_year" "$end_year" | { time /bin/bash $max_temp_preprocess_script > $seq_output ; } 2> >(tee "${seq_preprocess_time}" >&2)

## TODO: Figure a way to redirect stderr from curl not having input
echo "Executing the preprocessing max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_preprocess_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_preprocess_time}" >&2)
## This equivalence takes a very long time to check (uncomment with caution)
# echo "Checking for output equivalence..."
# diff -s $seq_output $pash_width_16_output | head

## Copy the sequential preprocess output to another file so that it doesn't get overwritten
echo "Copying intermediate file..."
split -n l/16 -d "$seq_output" ${preprocess_output}_16_

## Export the input variable for the process script
export IN="${preprocess_output}_16_*"

echo "Executing the processing max-temp script with bash..."
{ time /bin/bash $max_temp_process_script > $seq_output ; } 2> >(tee "${seq_process_time}" >&2)

## TODO: There is a bug in expansion that leads to quotes being added to $IN expansion
echo "Executing the processing max-temp script with pash -w 16 (log in: /tmp/pash_16_log)"
{ time $PASH_TOP/pa.sh -w 16 --log_file /tmp/pash_16_log --output_time $max_temp_process_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_process_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
