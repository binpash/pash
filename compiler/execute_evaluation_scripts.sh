#!/bin/bash

n_inputs=(
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
        ./execute_compile_evaluation_script.sh "${microbenchmark}_${n_in}"
    done
done
