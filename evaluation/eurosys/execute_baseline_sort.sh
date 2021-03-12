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
        h) echo "There are two possible execution levels:"
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
    intermediary_prefix=""
else
    echo "Unrecognizable execution level: $evaluation_level"
    exit 1
fi

eval_directory="$PASH_TOP/evaluation/"
intermediary_dir="$PASH_TOP/evaluation/intermediary/"
script_dir="${eval_directory}scripts/"
microbenchmarks_dir="$PASH_TOP/evaluation/microbenchmarks/"
results="${eval_directory}results/baseline_sort/"

mkdir -p $results

for n_in in "${n_inputs[@]}"; do
    experiment="baseline_sort_${n_in}"
    sort_parallel_script="${intermediary_dir}sort_${n_in}_seq.sh"
    env_file="${intermediary_dir}sort_${n_in}_env.sh"

    echo "Generating input and intermediary scripts... be patient..."
    python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
                $script_dir "sort" $n_in $intermediary_dir $env_suffix

    . $env_file
    export $(cut -d= -f1 $env_file)

    p_n_in="$(( $n_in * 2 ))"
    experiment="baseline_sort_${intermediary_prefix}${p_n_in}"
    echo "Executing sort with parallel flag for parallelism: ${p_n_in}"
    { time /bin/bash $sort_parallel_script "${p_n_in}" > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_parallel.time" >&2)

    echo "Generating input and intermediary scripts... be patient..."
    python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
                $microbenchmarks_dir "sort" $n_in $intermediary_dir $env_suffix

    exec_script="${intermediary_dir}sort_${n_in}_seq.sh"
    experiment="baseline_sort_${intermediary_prefix}${n_in}"

    if [ "$n_in" -eq 2 ]; then
        echo "Executing sort with bash"
        { time /bin/bash $exec_script ; } 1> /tmp/bash_output 2> >(tee "${results}${experiment}_seq.time" >&2)
    fi

    echo "Executing pash (no eager) on sort with --width ${n_in}"
    { time $PASH_TOP/pa.sh -w "${n_in}" --log_file /tmp/pash_log --output_time --no_eager $exec_script ; } 1> /tmp/pash_output 2> >(tee "${results}${experiment}_pash_no_eager.time" >&2)
    diff -s /tmp/seq_output /tmp/pash_output | head

    echo "Executing pash on sort with --width ${n_in}"
    { time $PASH_TOP/pa.sh -w "${n_in}" --log_file /tmp/pash_log --output_time $exec_script ; } 1> /tmp/pash_output 2> >(tee "${results}${experiment}_pash.time" >&2)
    diff -s /tmp/seq_output /tmp/pash_output | head
done
