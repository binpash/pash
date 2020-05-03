#!/bin/bash

eval_directory="../evaluation/"
intermediary_dir="../evaluation/intermediary/"
script_dir="${eval_directory}scripts/"
results="${eval_directory}results/baseline_sort/"

mkdir -p $results

n_inputs=(
    1
    2
    4
    8
    16
    32
    64
    96
    128
)

for n_in in "${n_inputs[@]}"; do
    experiment="baseline_sort_${n_in}"
    exec_script="${intermediary_dir}sort_${n_in}_seq.sh"
    env_file="${intermediary_dir}sort_${n_in}_env.sh"

    echo "Generating input and intermediary scripts... be patient..."
    python3 generate_microbenchmark_intermediary_scripts.py \
            $script_dir "sort" $n_in $intermediary_dir

    . $env_file
    export $(cut -d= -f1 $env_file)

    echo "Executing sort with parallel flag for parallelism: ${n_in}"
    { time /bin/bash $exec_script "${n_in}" > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)
done
