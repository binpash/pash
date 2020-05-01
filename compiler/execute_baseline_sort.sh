#!/bin/bash

eval_directory="../evaluation/"
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
    seq_script="${eval_directory}scripts/sort.sh"

    echo "Executing sort with parallel flag for parallelism: ${n_in}"
    { time /bin/bash $seq_script ; } 2> >(tee "${results}${experiment}_seq.time" >&2)
done
