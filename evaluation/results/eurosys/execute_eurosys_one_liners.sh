#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## This sets up to what extent we run the evaluation.
## There are 3 levels:
## 1. Small inputs | --width 2, 16 | Only full PaSh config
## 2. Small inputs | --width 2, 16 | All PaSh configs
## 3. Big inputs | -- width 2, 4, 8, 16, 32, 64 | All PaSh configs
##
## Note that for the small inputs there could be some variance with the results
## (especially with higher widths).
evaluation_level=1

while getopts 'smlh' opt; do
    case $opt in
        s) evaluation_level=1 ;;
        m) evaluation_level=2 ;;
        l) evaluation_level=3 ;;
        h) echo "There are three possible execution levels:"
           echo "option -s: Small inputs | --width 2, 16 | Only full PaSh config"
           echo "option -m: Small inputs | --width 2, 16 | All PaSh configs"
           echo "option -l: Big inputs | -- width 2, 4, 8, 16, 32, 64 | All PaSh configs"
           exit 0 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

## TODO: Add a script that runs the parallel sort evaluation

if [ "$evaluation_level" -eq 1 ]; then
    echo "Executing small evaluation..."
    n_inputs=(
        2
        16
    )
    result_subdir="eurosys_small"
    env_suffix="small"
    intermediary_prefix="small_"
    microbenchmarks=(
        'minimal_grep;-a'       # EuroSys: nfa-regex
        'minimal_sort;-a'       # EuroSys: sort
        'topn;-a'               # EuroSys: top-n
        'wf;-a'                 # EuroSys: wf
        'spell;-a'              # EuroSys: spell
        'diff;-a'               # EuroSys: difference
        'bigrams;-a'            # EuroSys: bi-grams
        'set-diff;-a'           # EuroSys: set-difference
        'double_sort;-a'        # EuroSys: sort-sort
        'shortest_scripts;-a'   # EuroSys: shortest-scripts
    )
elif [ "$evaluation_level" -eq 2 ]; then
    echo "Executing medium evaluation..."
    n_inputs=(
        2
        16
    )
    result_subdir="eurosys_small"
    env_suffix="small"
    intermediary_prefix="small_"
    microbenchmarks=(
        'minimal_grep;-n;-a'        # EuroSys: nfa-regex
        'minimal_sort;;-n;-a'       # EuroSys: sort
        'topn;;-n;-a'               # EuroSys: top-n
        'wf;;-n;-a'                 # EuroSys: wf
        'spell;-e;-a'               # EuroSys: spell
        'diff;;-n;-a'               # EuroSys: difference
        'bigrams;-e;-a'             # EuroSys: bi-grams
        'set-diff;;-n;-a'           # EuroSys: set-difference
        'double_sort;;-n;-e;-a'     # EuroSys: sort-sort
        'shortest_scripts;;-n;-a'   # EuroSys: shortest-scripts
    )
elif [ "$evaluation_level" -eq 3 ]; then
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
    result_subdir="eurosys_standard"
    env_suffix=""
    intermediary_prefix=""
    microbenchmarks=(
        'minimal_grep;-n;-a'        # EuroSys: nfa-regex
        'minimal_sort;;-n;-a'       # EuroSys: sort
        'topn;;-n;-a'               # EuroSys: top-n
        'wf;;-n;-a'                 # EuroSys: wf
        'spell;-e;-a'               # EuroSys: spell
        'diff;;-n;-a'               # EuroSys: difference
        'bigrams;-e;-a'             # EuroSys: bi-grams
        'set-diff;;-n;-a'           # EuroSys: set-difference
        'double_sort;;-n;-e;-a'     # EuroSys: sort-sort
        'shortest_scripts;;-n;-a'   # EuroSys: shortest-scripts
    )
else
    echo "Unrecognizable execution level: $evaluation_level"
    exit 1
fi

microbenchmarks_dir="$PASH_TOP/evaluation/microbenchmarks/"
intermediary_dir="$PASH_TOP/evaluation/${intermediary_prefix}intermediary/"
mkdir -p $intermediary_dir
mkdir -p "$PASH_TOP/evaluation/results/$result_subdir/"

for microbenchmark_config in "${microbenchmarks[@]}"; do
    IFS=";" read -r -a flags <<< "${microbenchmark_config}"
    microbenchmark=${flags[0]}
    echo "Executing: $microbenchmark"
    # Execute the sequential script on the first run only
    exec_seq="-s"
    for n_in in "${n_inputs[@]}"; do

        ## Generate the intermediary script
        echo "Generating input and intermediary scripts... be patient..."
        python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
                $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir $env_suffix

        for flag in "${flags[@]:1}"; do
            echo "Flag: ${flag}"

            ## Execute the intermediary script
            "$PASH_TOP/evaluation/execute_compile_evaluation_script.sh" $exec_seq $flag "${microbenchmark}" "${n_in}" $result_subdir $intermediary_prefix > /dev/null 2>&1

            ## Only run the sequential the first time around
            exec_seq=""
        done
    done
done
