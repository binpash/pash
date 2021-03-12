#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export DEBUG=0
# export DEBUG=1 # Uncomment to print pash output

microbenchmarks_dir="${PASH_TOP}/evaluation/microbenchmarks/"
intermediary_dir="${PASH_TOP}/evaluation/test_intermediary/"
test_results_dir="${PASH_TOP}/evaluation/results/test_results/"

echo "Deleting eager intermediate files..."
rm -f /tmp/eager*
mkdir -p $intermediary_dir
rm -rf "$test_results_dir"
mkdir -p "$test_results_dir"

n_inputs=(
    2
    8
)

configurations=(
    ""
    "--r_split"
)

## Tests where the compiler will not always succeed (e.g. because they have mkfifo)
script_microbenchmarks=(
    diff                 # (quick-abort) BUG: Might have to do with the named pipes, and the fact that they are reused for parallel and sequential script.
    set-diff             # TODO: Handle redirection after reduce
    export_var_script    # Tests whether exported variables in the scripts that are processed by PaSh runtime are visible to the rest of the script.
    comm-par-test        # Small comm test to ensure non-parallelizability
    comm-par-test2       # Small comm test with input redirection and hyphen
    tee_web_index_bug    # Tests a tee bug from web index
)

pipeline_microbenchmarks=(
    grep                 # One-liner
    minimal_sort         # One-liner
    minimal_grep         # One-liner
    topn                 # One-liner
    wf                   # One-liner
    spell                # One-liner
    shortest_scripts     # One-liner
    bigrams              # One-liner
    alt_bigrams          # One-liner
    deadlock_test        # Test to check deadlock prevention using drain_stream
    double_sort          # Checks maximum peformance gains from split
    no_in_script         # Tests whether a script can be executed by our infrastructure without having its input in a file called $IN
    for_loop_simple      # Tests whether PaSh can handle a for loop where the body is parallelizable
    minimal_grep_stdin   # Tests whether PaSh can handle a script that reads from stdin
    # micro_1000           # Not being run anymore, as it is very slow. Tests whether the compiler is fast enough. It is a huge pipeline without any computation.
    micro_10           # A small version of the pipeline above for debugging.
    sed-test             # Tests all sed occurences in our evaluation to make sure that they work
    fun-def              # Tests whether PaSh can handle a simple function definition
    tr-test              # Tests all possible behaviors of tr that exist in our evaluation
)

execute_pash_and_check_diff() {
    if [ "$DEBUG" -eq 1 ]; then
        { time "$PASH_TOP/pa.sh" $@ ; } 1> "$pash_output" 2> >(tee "${pash_time}" >&2) &&
        diff -s "$seq_output" "$pash_output" | head | tee -a "${pash_time}" >&2
    else
        { time "$PASH_TOP/pa.sh" $@ ; } 1> "$pash_output" 2> "${pash_time}" &&
        diff -s "$seq_output" "$pash_output" | head >> "${pash_time}"
    fi
}

execute_tests() {
    assert_correctness="$1"
    microbenchmarks=("${@:2}")

    microbenchmark_configs=( )
    for i in "${!microbenchmarks[@]}"; do
        all_flags=${test_flags[@]}
        microbenchmark_configs[$i]="${microbenchmarks[$i]};${all_flags// /;}"
    done

    ## This is almost the same loop as the one in execute_evaluation_scripts
    for microbenchmark_config in "${microbenchmark_configs[@]}"; do
        IFS=";" read -r -a flags <<< "${microbenchmark_config}"
        microbenchmark=${flags[0]}
        echo "Executing test: $microbenchmark"
        # Execute the sequential script on the first run only
        
        prefix="${microbenchmarks_dir}/${microbenchmark}"

        export seq_output="${intermediary_dir}/${microbenchmark}_seq_output"
        seq_time="$test_results_dir/${microbenchmark}_seq.time"

        script_to_execute="${prefix}.sh"
        env_file="${prefix}_env_test.sh"
        funs_file="${prefix}_funs.sh"
        input_file="${prefix}_test.in"

        if [ -f "$env_file" ]; then
            . $env_file
            vars_to_export=$(cut -d= -f1 $env_file)
            if [ ! -z "$vars_to_export" ]; then
                export $vars_to_export
            fi
        else
            echo "|-- Doesn't have env file"
        fi

        ## Export necessary functions
        if [ -f "$funs_file" ]; then
            source $funs_file
        fi

        ## Redirect the input if there is an input file
        stdin_redir="/dev/null"
        if [ -f "$input_file" ]; then
            stdin_redir="$(cat "$input_file")"
            echo "|-- Has input file: $stdin_redir"
        fi

        echo "|-- Executing the script with bash..."
        cat $stdin_redir | { time /bin/bash "$script_to_execute" > $seq_output ; } 2> "${seq_time}"

        for conf in "${configurations[@]}"; do
            for n_in in "${n_inputs[@]}"; do
                echo "|-- Executing with pash --width ${n_in} ${conf}..."
                export pash_time="${test_results_dir}/${microbenchmark}_${n_in}_distr_${conf}.time"
                export pash_output="${intermediary_dir}/${microbenchmark}_${n_in}_pash_output"

                cat $stdin_redir |
                    execute_pash_and_check_diff -d 1 $assert_correctness ${conf} --width "${n_in}" --output_time $script_to_execute
            done
        done
    done
}

execute_tests "" "${script_microbenchmarks[@]}"
execute_tests "--assert_compiler_success" "${pipeline_microbenchmarks[@]}"

echo "Below follow the identical outputs:"
grep --files-with-match "are identical" "$test_results_dir"/*_distr*.time |
    sed "s,^$PASH_TOP/,,"

echo "Below follow the non-identical outputs:"
grep -L "are identical" "$test_results_dir"/*_distr*.time |
    sed "s,^$PASH_TOP/,,"

TOTAL_TESTS=$(ls -la "$test_results_dir"/*_distr*.time | wc -l)
PASSED_TESTS=$(grep --files-with-match "are identical" "$test_results_dir"/*_distr*.time | wc -l)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
