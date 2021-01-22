#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## This sets up to what extent we run the evaluation.
## There are 2 levels:
## 1. Small input | --width 2, 16
## 2. Big input | -- width 2, 4, 8, 16, 32, 64
evaluation_level=1

while getopts 'slh' opt; do
    case $opt in
        s) evaluation_level=1 ;;
        l) evaluation_level=2 ;;
        h) echo "There are three possible execution levels:"
           echo "option -s: Small input | --width 2, 16"
           echo "option -l: Big input | -- width 2, 4, 8, 16, 32, 64"
           exit 0 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

if [ "$evaluation_level" -eq 1 ]; then
    echo "Executing small baseline sort evaluation..."
    n_inputs=(
        2
        16
    )
    env_suffix="small"
    env_suffix_underscore="_small"
    intermediary_prefix="small_"
elif [ "$evaluation_level" -eq 2 ]; then
    echo "Executing large baseline sort evaluation..."
    n_inputs=(
        2
        4
        8
        16
        32
        64
    )
    env_suffix=""
    env_suffix_underscore=""
    intermediary_prefix=""
else
    echo "Unrecognizable execution level: $evaluation_level"
    exit 1
fi

eval_directory="$PASH_TOP/evaluation/"
intermediary_dir="$PASH_TOP/evaluation/intermediary/"
script_dir="${eval_directory}scripts/"
results="${eval_directory}results/baseline_sort/"

mkdir -p $results

for n_in in "${n_inputs[@]}"; do
    experiment="baseline_sort_${n_in}"
    echo "$experiment"
    exec_script="${intermediary_dir}sort_${n_in}_seq.sh"
    env_file="${intermediary_dir}sort_${n_in}_env${env_suffix_underscore}.sh"

    echo "Generating input and intermediary scripts... be patient..."
    python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
                $script_dir "sort" $n_in $intermediary_dir $env_suffix

    . $env_file
    export $(cut -d= -f1 $env_file)

    echo "Executing sort with parallel flag for parallelism: ${n_in}"
    { time /bin/bash $exec_script "${n_in}" > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_${intermediary_prefix}seq.time" >&2)

    echo "Executing pash on sort with --width ${n_in}"
    { time $PASH_TOP/pa.sh -w "${n_in}" --log_file /tmp/pash_2_log --output_time $exec_script ; } 1> /tmp/pash_output 2> >(tee "${results}${experiment}_${intermediary_prefix}pash.time" >&2)
done
