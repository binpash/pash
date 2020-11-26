#!/bin/bash

microbenchmarks_dir="../evaluation/microbenchmarks/"
intermediary_dir="../evaluation/intermediary/"

echo "Deleting eager intermediate files..."
rm -f /tmp/eager*
mkdir -p $intermediary_dir

gnu_parallel_flag=0

n_inputs=(
    2
    4
    8
    16
    32
    64
)

microbenchmarks=(
    'grep;-n;-e'           # PLDI
    'minimal_sort;;-n;-e'         # PLDI
    'minimal_grep;-n;-e'         # PLDI
    'topn;;-n;-e'                 # PLDI
    'wf;;-n;-e'                   # PLDI
    'spell;-e;-a'                # PLDI
    'shortest_scripts;;-n;-e'     # PLDI
    # micro_1000           # PLDI

    'bigrams;-e;-a'              # TODO: Fix bug. Run with good split.
    'alt_bigrams;;-n;-e'          # Optimized version of Bigrams
    'diff;;-n;-e'                 # TODO: Optimize diff
    'set-diff;;-n;-e'             # TODO: Handle redirection after reduce

    # sort                 # For comparison with sort --parallel

    ## Tests
    # deadlock_test          # Test to check deadlock prevention using drain_stream
    'double_sort;;-n;-e;-a'            # Checks maximum peformance gains from split

    ## TODO: Add some more for OSDI
    # wc                   # Extremely simple
    # page-count           # TODO: Change it so that it is one pipeline and it has many files
    # sieve                # Nice to show posix compliance
    # genquality           # TODO: Nice big pipeline. Handle this
    # genome-diff          # TODO: Nice big pipeline. Handle this
    # symtab               # Probably no benefit (since it runs some sequential pure function)
    # for_loop_simple      # Probably not for OSDI
)

for microbenchmark_config in "${microbenchmarks[@]}"; do
    IFS=";" read -r -a flags <<< "${microbenchmark_config}"
    microbenchmark=${flags[0]}
    echo "Executing: $microbenchmark"
    # Execute the sequential script on the first run only
    exec_seq="-s"
    for n_in in "${n_inputs[@]}"; do

        ## Generate the intermediary script
        echo "Generating input and intermediary scripts... be patient..."
        python3 generate_microbenchmark_intermediary_scripts.py \
                $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir

        for flag in "${flags[@]:1}"; do
            echo "Flag: ${flag}"

            ## Execute the intermediary script
            ./execute_compile_evaluation_script.sh $exec_seq $flag "${microbenchmark}" "${n_in}"
            rm -f /tmp/eager*

            ## Only run the sequential the first time around
            exec_seq=""
        done

        ## Execute the gnu parallel
        if [ "$gnu_parallel_flag" -eq 1 ]; then
            ./execute_gnu_parallel_script.sh $microbenchmark $n_in
        fi
    done
done
