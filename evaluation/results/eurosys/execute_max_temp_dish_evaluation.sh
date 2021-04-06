#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## There are two possible execution levels:
## option -s: end_year=2000
## option -l: end_year=2004 (The EuroSys evaluation)
## option -e: Run extended (separate preprocessing from processing) (The EuroSys evaluation)
start_year=2000
end_year=2000 # For the small evaluation
execute_separate_flag=0 # Whether to execute processing and preprocessing separately

while getopts 'sleh' opt; do
    case $opt in
        s) end_year=2000 ;;
        l) end_year=2004 ;;
        e) execute_separate_flag=1 ;;
        h) echo "There are three possible execution levels:"
           echo "option -s: end_year=2000"
           echo "option -l: end_year=2004 (The EuroSys evaluation)"
           echo "option -e: Run extended (separate preprocessing from processing) (The EuroSys evaluation)"
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

temp_dir=max_temp_tmp_results
mkdir -p $temp_dir

seq_output="${temp_dir}/max_temp_seq_output"
pash_width_16_output="${temp_dir}/max_temp_pash_16_output"
seq_time="$results_dir/max-temp-complete-$start_year-$end_year-seq.time"
pash_width_16_time="$results_dir/max-temp-complete-$start_year-$end_year-16-pash.time"

echo "Executing the complete max-temp script with bash..."
seq "$start_year" "$end_year" | { time /bin/bash $max_temp_complete_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

echo "Executing the complete max-temp script with pash -w 16 (log in: ${temp_dir}/max_temp_pash_16_log)"
seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file "${temp_dir}/max_temp_pash_16_log" --output_time $max_temp_complete_script ; } 1> "$pash_width_16_output" 2> >(tee "$pash_width_16_time" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
 
if [ "$execute_separate_flag" -eq 1 ]; then
    echo "Extended: Executing preprocessing and processing separately"

    max_temp_preprocess_script="${eval_dir}/scripts/max-temp-preprocess.sh"
    max_temp_process_script="${eval_dir}/scripts/max-temp-process.sh"

    seq_preprocess_time="$results_dir/max-temp-preprocess-$start_year-$end_year-seq.time"
    pash_width_16_preprocess_time="$results_dir/max-temp-preprocess-$start_year-$end_year-16-pash.time"
    seq_process_time="$results_dir/max-temp-process-$start_year-$end_year-seq.time"
    pash_width_16_process_time="$results_dir/max-temp-process-$start_year-$end_year-16-pash.time"
    preprocess_output="${temp_dir}/max-temp-preprocess-output"

    echo "Executing the preprocessing max-temp script with bash..."
    seq "$start_year" "$end_year" | { time /bin/bash $max_temp_preprocess_script > $seq_output ; } 2> >(tee "${seq_preprocess_time}" >&2)

    echo "Executing the preprocessing max-temp script with pash -w 16 (log in: ${temp_dir}/pash_16_log)"
    seq "$start_year" "$end_year" | { time $PASH_TOP/pa.sh -w 16 --log_file "${temp_dir}/pash_16_log" --output_time $max_temp_preprocess_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_preprocess_time}" >&2)
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

    echo "Executing the processing max-temp script with pash -w 16 (log in: ${temp_dir}/pash_16_log)"
    { time $PASH_TOP/pa.sh -w 16 --log_file "${temp_dir}/pash_16_log" --output_time $max_temp_process_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_process_time}" >&2)
    echo "Checking for output equivalence..."
    diff -s $seq_output $pash_width_16_output | head
fi
