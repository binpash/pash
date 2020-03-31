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
    grep
    minimal_sort
    minimal_grep         # PLDI
    topn                 # PLDI
    wf                   # PLDI
    spell                # PLDI
    shortest_scripts     # PLDI
    micro_1000           # PLDI

    # bigrams              # TODO: Fix bug. Run with good split.
    alt_bigrams          # Not so true (as many steps are combined for better MapReduce)

    ## TODO: Add some more for OSDI
    # wc                   # Extremely simple
    diff                 # TODO: Optimize diff
    # set-diff             # TODO: Handle that
    # page-count           # TODO: Change it so that it is one pipeline and it has many files
    # sieve                # Nice to show posix compliance
    # genquality           # TODO: Nice big pipeline. Handle this
    # genome-diff          # TODO: Nice big pipeline. Handle this
    # symtab               # Probably no benefit (since it runs some sequential pure function)
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
