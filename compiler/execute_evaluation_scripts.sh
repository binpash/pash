#!/bin/bash

microbenchmarks_dir="../evaluation/microbenchmarks/"
intermediary_dir="../evaluation/intermediary/"

echo "Deleting eager intermediate files..."
rm -f /tmp/eager*
mkdir -p $intermediary_dir

n_inputs=(
    # 1
    2
    4
    8
    16
    32
    64
    96
    128
    # 10
    # 20
    # 50
    # 100
)

microbenchmarks=(
    grep                 # PLDI
    minimal_sort         # PLDI
    minimal_grep         # PLDI
    topn                 # PLDI
    wf                   # PLDI
    spell                # PLDI
    shortest_scripts     # PLDI
    # micro_1000           # PLDI

    # bigrams              # TODO: Fix bug. Run with good split.
    alt_bigrams          # Optimized version of Bigrams
    diff                 # TODO: Optimize diff
    set-diff             # TODO: Handle redirection after reduce

    ## Tests
    # deadlock_test          # Test to check deadlock prevention using drain_stream

    ## TODO: Add some more for OSDI
    # wc                   # Extremely simple
    # page-count           # TODO: Change it so that it is one pipeline and it has many files
    # sieve                # Nice to show posix compliance
    # genquality           # TODO: Nice big pipeline. Handle this
    # genome-diff          # TODO: Nice big pipeline. Handle this
    # symtab               # Probably no benefit (since it runs some sequential pure function)
    # for_loop_simple      # Probably not for OSDI
)

## Note: Maybe we need to tune the configuration (fan-out, batch-size)
##       before some specific benchmarks
for microbenchmark in "${microbenchmarks[@]}"; do
    # Execute the sequential script on the first run only
    exec_seq="-s"
    for n_in in "${n_inputs[@]}"; do

        ## Generate the intermediary script
        echo "Generating input and intermediary scripts... be patient..."
        python3 generate_microbenchmark_intermediary_scripts.py \
                $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir

        ## Execute the intermediary script with eager
        ./execute_compile_evaluation_script.sh $exec_seq -e "${microbenchmark}" "${n_in}"
        rm -f /tmp/eager*

        # Only execute the sequential once
        exec_seq=""

        ## Execute the intermediary script without eager
        ./execute_compile_evaluation_script.sh $exec_seq "${microbenchmark}" "${n_in}"

        ## Execute the intermediary script with the naive eager
        ./execute_compile_evaluation_script.sh $exec_seq -n "${microbenchmark}" "${n_in}"
        rm -f /tmp/eager*
    done
done
