#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## There are two possible execution levels:
## options -s: 1,000 urls (about 1.5 minutes in bash)
## options -l: 100,000 urls (a couple hours in bash)
input_number=1000 # About 1.5 minutes in bash
# input_number=100 # About 7 seconds in bash

while getopts 'slh' opt; do
    case $opt in
        s) input_number=1000 ;;
        l) input_number=100000 ;;
        h) echo "There are two possible execution levels:"
           echo "option -s: 1,000 urls (about 1.5 minutes in bash)"
           echo "option -l: 100,000 urls (a couple hours in bash) (EuroSys evaluation)"
           exit 0 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

eval_dir="$PASH_TOP/evaluation/"
directory="${eval_dir}/scripts/web-index/"
results_dir="${eval_dir}/results/"
input_dir="${HOME}/wikipedia/"

export IN="$input_dir/index_h_${input_number}.txt"
export WIKI="${input_dir}"
export WEB_INDEX_DIR="${directory}"

web_index_script="${eval_dir}/scripts/web-index.sh"

temp_dir=web_index_tmp_results
mkdir -p "$temp_dir"

seq_output="${temp_dir}/seq_output"
pash_width_2_output="${temp_dir}/pash_2_output"
pash_width_16_output="${temp_dir}/pash_16_output"
seq_time="$results_dir/web-index-${input_number}-seq.time"
pash_width_2_time="$results_dir/web-index-${input_number}-2-pash.time"
pash_width_16_time="$results_dir/web-index-${input_number}-16-pash.time"

echo "Executing the script with bash..."
{ time /bin/bash $web_index_script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

echo "Executing the script with pash -w 2 (log in: ${temp_dir}/pash_2_log)"
{ time $PASH_TOP/pa.sh -w 2 --log_file "${temp_dir}/pash_2_log" --output_time $web_index_script ; } 1> "$pash_width_2_output" 2> >(tee "${pash_width_2_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_2_output | head

echo "Executing the script with pash -w 16 (log in: ${temp_dir}/pash_16_log)"
{ time $PASH_TOP/pa.sh -w 16 --log_file "${temp_dir}/pash_16_log" --output_time $web_index_script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_time}" >&2)
echo "Checking for output equivalence..."
diff -s $seq_output $pash_width_16_output | head
