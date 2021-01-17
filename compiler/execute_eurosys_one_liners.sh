#!/bin/bash

## TODO: Set this up with an input argument
small_evaluation_flag=1

## TODO: Make sure that the pash results in the plotting script are collected from distr-auto-split.time

## TODO: Add a script that runs the parallel sort evaluation

if [ "$small_evaluation_flag" -eq 1 ]; then
    echo "Executing small evaluation..."
    n_inputs=(
        2
        16
    )
    result_subdir="eurosys_small"
    env_suffix="small"
    intermediary_prefix="small_"
else
    echo "Executing standard evaluation..."
    n_inputs=(
        2
        4
        8
        16
        32
        64
    )
    ## TODO: Maybe change the result_subdir for the full evaluation
    result_subdir=""
    env_suffix=""
    intermediary_prefix=""
fi

microbenchmarks=(
    'minimal_grep;-n;-e'        # EuroSys: nfa-regex
    'minimal_sort;;-n;-e'       # EuroSys: sort
    'topn;;-n;-e'               # EuroSys: top-n
    'wf;;-n;-e'                 # EuroSys: wf
    'spell;-e;-a'               # EuroSys: spell
    'diff;;-n;-e'               # EuroSys: difference
    'bigrams;-e;-a'             # EuroSys: bi-grams
    'set-diff;;-n;-e'           # EuroSys: set-difference
    'double_sort;;-n;-e;-a'     # EuroSys: sort-sort
    'shortest_scripts;;-n;-e'   # EuroSys: shortest-scripts
)

microbenchmarks_dir="../evaluation/microbenchmarks/"
intermediary_dir="../evaluation/${intermediary_prefix}intermediary/"
mkdir -p $intermediary_dir
mkdir -p "../evaluation/results/$result_subdir/"


echo "Deleting eager intermediate files..."
rm -f /tmp/eager*
mkdir -p $intermediary_dir

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
                $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir $env_suffix

        for flag in "${flags[@]:1}"; do
            echo "Flag: ${flag}"

            ## Execute the intermediary script
            ./execute_compile_evaluation_script.sh $exec_seq $flag "${microbenchmark}" "${n_in}" $result_subdir $intermediary_prefix
            rm -f /tmp/eager*

            ## Only run the sequential the first time around
            exec_seq=""
        done
    done
done
