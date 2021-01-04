#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
evaluation_dir="${PASH_TOP?Variable PASH_TOP is undefined}/evaluation/"
microbenchmarks_dir="${evaluation_dir}/microbenchmarks/"
results="${evaluation_dir}results/perf/"
mkdir -p "$results"

seq_output="/tmp/pash_evaluation_seq_output"
pash_output="/tmp/pash_evaluation_output"

bash="/bin/bash"
pash="$PASH_TOP/pa.sh --output_time"

echo "Deleting eager intermediate files..."
rm -f /tmp/eager*

widths=(
    2
)

microbenchmarks=(
    # 'grep'           # PLDI
    # 'minimal_sort'         # PLDI
    # 'minimal_grep'         # PLDI
    # 'topn'                 # PLDI
    # 'wf'                   # PLDI
    # 'spell'                # PLDI
    # 'shortest_scripts'     # PLDI
    # # micro_1000           # PLDI

    # 'bigrams'              # TODO: Fix bug. Run with good split.
    # 'alt_bigrams'          # Optimized version of Bigrams
    # 'diff'                 # TODO: Optimize diff
    # 'set-diff'             # TODO: Handle redirection after reduce

    # ## Tests
    'double_sort'            # Checks maximum peformance gains from split

    ## TODO: Add some more for OSDI
    # wc                   # Extremely simple
    # page-count           # TODO: Change it so that it is one pipeline and it has many files
    # sieve                # Nice to show posix compliance
    # genquality           # TODO: Nice big pipeline. Handle this
    # genome-diff          # TODO: Nice big pipeline. Handle this
    # symtab               # Probably no benefit (since it runs some sequential pure function)
    # for_loop_simple      # Probably not for OSDI
)

##
## Scripts to execute microbenchmarks
##

exec_microbenchmark()
{
    shell="$1"
    experiment="$2"
    output="$3"
    result_file="$4"
    prefix="${microbenchmarks_dir}${experiment}"
    env_file="${prefix}_env.sh"
    funs_file="${prefix}_funs.sh"
    script="${prefix}.sh"
    input_file="${prefix}.in"

    . $env_file
    vars_to_export=$(cut -d= -f1 $env_file)
    if [ ! -z "$vars_to_export" ]; then
        export $vars_to_export
    fi

    ## Export necessary functions
    if [ -f "$funs_file" ]; then
        source $funs_file
    fi

    ## Redirect the input if there is an input file
    stdin_redir="/dev/null"
    if [ -f "$input_file" ]; then
        stdin_redir="$(cat "$input_file")"
    fi

    cat $stdin_redir | { time $shell $script ; } 1> $output 2> >(tee "$result_file" >&2)
}

exec_seq_microbenchmark()
{
    experiment="$1"
    result_file="${results}${experiment}_seq.time"
    echo "Executing sequential microbenchmark: $experiment "
    exec_microbenchmark "$bash" "$experiment" "$seq_output" "$result_file"
}

exec_pash_microbenchmark()
{
    experiment="$1"
    width="$2"
    result_file="${results}${experiment}_pash.time"
    echo "Executing PaSh microbenchmark: $experiment "
    exec_microbenchmark "$pash -w $width" "$experiment" "$pash_output" "$result_file"
    echo "Checking for equivalence..." &&
    diff -s "$seq_output" "$pash_output" | head | tee -a "$result_file"
}

##
## This is the main loop that executes all programs with all configurations
##

for microbenchmark in "${microbenchmarks[@]}"; do
    echo "Executing: $microbenchmark"

    exec_seq_microbenchmark "$microbenchmark"
    for width in "${widths[@]}"; do
        exec_pash_microbenchmark "$microbenchmark" "$width"
    done
done
