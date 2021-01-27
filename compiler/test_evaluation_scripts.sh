#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

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
    # 8
)

## Tests where the compiler will not always succeed (e.g. because they have mkfifo)
script_microbenchmarks=(
    # diff                 # (quick-abort) BUG: Might have to do with the named pipes, and the fact that they are reused for parallel and sequential script.
    # set-diff             # TODO: Handle redirection after reduce
    # export_var_script    # Tests whether exported variables in the scripts that are processed by PaSh runtime are visible to the rest of the script.
    # comm-par-test        # Small comm test to ensure non-parallelizability
    # comm-par-test2       # Small comm test with input redirection and hyphen
    # tee_web_index_bug    # Tests a tee bug from web index
)

pipeline_microbenchmarks=(
    # grep                 # One-liner
    minimal_sort         # One-liner
    # minimal_grep         # One-liner
    # topn                 # One-liner
    # wf                   # One-liner
    # spell                # One-liner
    # shortest_scripts     # One-liner
    # bigrams              # One-liner
    # alt_bigrams          # One-liner
    # deadlock_test        # Test to check deadlock prevention using drain_stream
    # double_sort          # Checks maximum peformance gains from split
    # no_in_script         # Tests whether a script can be executed by our infrastructure without having its input in a file called $IN
    # for_loop_simple      # Tests whether PaSh can handle a for loop where the body is parallelizable
    # minimal_grep_stdin   # Tests whether PaSh can handle a script that reads from stdin
    # micro_1000           # Tests whether the compiler is fast enough. It is a huge pipeline without any computation.
    # sed-test             # Tests all sed occurences in our evaluation to make sure that they work
    # fun-def              # Tests whether PaSh can handle a simple function definition
    # tr-test              # Tests all possible behaviors of tr that exist in our evaluation
)

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

        
        prefix="${intermediary_dir}/${microbenchmark}"

        seq_output="${intermediary_dir}/${microbenchmark}_seq_output"
        # pash_width_2_output=/tmp/pash_2_output
        seq_time="$test_results_dir/${microbenchmark}_seq.time"
        # pash_width_2_time="$results_dir/web-index-${input_number}-2-pash.time"

        # echo "Executing the script with pash -w 2 (log in: /tmp/pash_2_log)"
        # { time $PASH_TOP/pa.sh -w 2 --log_file /tmp/pash_2_log --output_time $web_index_script ; } 1> "$pash_width_2_output" 2> >(tee "${pash_width_2_time}" >&2)
        # echo "Checking for output equivalence..."
        # diff -s $seq_output $pash_width_2_output | head




        exec_seq="-s"
        for n_in in "${n_inputs[@]}"; do
            echo "Number of inputs: ${n_in}"

            ## TODO: Similarly stop using this...
            ## Generate the intermediary script
            python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
                    $microbenchmarks_dir $microbenchmark $n_in $intermediary_dir "test"

            ## TODO: Move these out once we remove the generate...
            script_to_execute="${prefix}_${n_in}_seq.sh"
            env_file="${prefix}_${n_in}_env.sh"
            funs_file="${prefix}_${n_in}_funs.sh"
            input_file="${prefix}_${n_in}.in"
            echo "$script_to_execute"

            if [ "$exec_seq" == "-s" ]; then
                echo "Environment:"
                # cat $env_file
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
                    echo "Has input file: $stdin_redir"
                fi


                echo "Executing the script with bash..."
                cat $stdin_redir | { time /bin/bash "$script_to_execute" > $seq_output ; } 2> >(tee "${seq_time}" >&2)
            fi
            ## TODO: Stop using execute_compile_evaluation_script.sh

            ## Execute the intermediary script
            "$PASH_TOP/evaluation/execute_compile_evaluation_script.sh" $assert_correctness -a "${microbenchmark}" "${n_in}" "test_results" "test_" # > /dev/null 2>&1
            rm -f /tmp/eager*

            ## Only run the sequential the first time around
            exec_seq=""
        done
    done
}

execute_tests "" "${script_microbenchmarks[@]}"
execute_tests "-c" "${pipeline_microbenchmarks[@]}"

echo "Below follow the identical outputs:"
grep --files-with-match "are identical" "$test_results_dir"/*_distr*.time

echo "Below follow the non-identical outputs:"
grep -L "are identical" "$test_results_dir"/*_distr*.time

TOTAL_TESTS=$(ls -la "$test_results_dir"/*_distr*.time | wc -l)
PASSED_TESTS=$(grep --files-with-match "are identical" "$test_results_dir"/*_distr*.time | wc -l)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
