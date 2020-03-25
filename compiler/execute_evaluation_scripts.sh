#!/bin/bash

microbenchmarks_dir="../evaluation/microbenchmarks/"
intermediary_dir="../evaluation/intermediary/"

mkdir -p $intermediary_dir

n_inputs=(
    1
    2
    4
    10
    20
    50
    100
)

microbenchmarks=(
    minimal_sort
    minimal_grep
    topn
    wf
    grep
    spell
    shortest_scripts
    micro_1000
)

## Note: Maybe we need to tune the configuration (fan-out, batch-size)
##       before some specific benchmarks
for microbenchmark in "${microbenchmarks[@]}"; do
    for n_in in "${n_inputs[@]}"; do

        ## Generate the intermediary script
        python3 generate_microbenchmark_intermediary_scripts.py \
                $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir

        ## Execute the intermediary script
        ./execute_compile_evaluation_script.sh "${microbenchmark}_${n_in}"
    done
done
