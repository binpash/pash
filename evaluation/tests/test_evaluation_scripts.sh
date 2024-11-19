#!/bin/bash
# time: print real in seconds, to simplify parsing
## Necessary to set PASH_TOP
cd $(dirname $0)
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export DEBUG=0
export PASH_LOG=1
# export DEBUG=1 # Uncomment to print pash output
## Determines whether the experimental pash flags will be tested. 
## By default they are not.
export EXPERIMENTAL=0
for item in $@
do
    if [ "--debug" == "$item" ] || [ "-d" == "$item" ]; then
        export DEBUG=1
    fi
    if [ "--no-pash-log" == "$item" ]; then
        export PASH_LOG=0
    fi
    if [ "--experimental" == "$item" ]; then
        export EXPERIMENTAL=1
    fi
done

microbenchmarks_dir="${PASH_TOP}/evaluation/tests"
intermediary_dir="${PASH_TOP}/evaluation/tests/test_intermediary"
test_results_dir="${PASH_TOP}/evaluation/tests/results"
results_time="$test_results_dir/results.time"
results_time_bash=${results_time}_bash
results_time_pash=${results_time}_pash

echo "Deleting eager intermediate files..."
rm -rf "$test_results_dir"
rm -rf "$intermediary_dir"
mkdir -p $intermediary_dir
mkdir -p "$test_results_dir"

echo "Generating inputs..."
cd "$microbenchmarks_dir/input"
./setup.sh
cd -

n_inputs=(
    2
    8
)

if [ "$EXPERIMENTAL" -eq 1 ]; then
    configurations=(
        ""
    )
else
    configurations=(
        "--profile_driven"
    )
fi


## Tests where the compiler will not always succeed (e.g. because they have mkfifo)
script_microbenchmarks=(
    diff                 # (quick-abort) BUG: Might have to do with the named pipes, and the fact that they are reused for parallel and sequential script.
    set-diff             # TODO: Handle redirection after reduce
    export_var_script    # Tests whether exported variables in the scripts that are processed by PaSh runtime are visible to the rest of the script.
    comm-par-test        # Small comm test to ensure non-parallelizability
    comm-par-test2       # Small comm test with input redirection and hyphen
    tee_web_index_bug    # Tests a tee bug from web index
    fun-def              # Tests whether PaSh can handle a simple function definition
    bigrams              # One-liner
    spell-grep           # Spell variant with `grep -f` instead of `comm`
)

pipeline_microbenchmarks=(
    grep                 # One-liner
    minimal_sort         # One-liner
    minimal_grep         # One-liner
    topn                 # One-liner
    wf                   # One-liner
    spell                # One-liner
    shortest_scripts     # One-liner
    alt_bigrams          # One-liner
    deadlock_test        # Test to check deadlock prevention using drain_stream
    double_sort          # Checks maximum peformance gains from split
    no_in_script         # Tests whether a script can be executed by our infrastructure without having its input in a file called $IN
    for_loop_simple      # Tests whether PaSh can handle a for loop where the body is parallelizable
    minimal_grep_stdin   # Tests whether PaSh can handle a script that reads from stdin
    micro_10             # A small version of the pipeline above for debugging.
    sed-test             # Tests all sed occurences in our evaluation to make sure that they work
    tr-test              # Tests all possible behaviors of tr that exist in our evaluation
    grep-test            # Tests some interesting grep invocations
    ann-agg              # Tests custom aggregators in annotations
    # # # # micro_1000           # Not being run anymore, as it is very slow. Tests whether the compiler is fast enough. It is a huge pipeline without any computation.
)



execute_pash_and_check_diff() {
    TIMEFORMAT="%3R" # %3U %3S"
    if [ "$DEBUG" -eq 1 ]; then
        { time "$PASH_TOP/pa.sh" $@ ; } 1> "$pash_output" 2> >(tee -a "${pash_time}" >&2) &&
        diff -s "$seq_output" "$pash_output" | head | tee -a "${pash_time}" >&2
    else

        { time "$PASH_TOP/pa.sh" $@ ; } 1> "$pash_output" 2>> "${pash_time}" &&
        b=$(cat "$pash_time"); 
        test_diff_ec=$(cmp -s "$seq_output" "$pash_output" && echo 0 || echo 1)
        # differ
        script=$(basename $script_to_execute)
        if [ $test_diff_ec -ne 0 ]; then
            c=$(diff -s "$seq_output" "$pash_output" | head)
            echo "$c$b" > "${pash_time}"
            echo "$script are not identical" >> $test_results_dir/result_status
        else
            echo "Files $seq_output and $pash_output are identical" > "${pash_time}"
            echo "$script are identical" >> $test_results_dir/result_status
        fi

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

        export script_to_execute="${prefix}.sh"
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
            echo "|-- Does not have env file"
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

        TIMEFORMAT="${microbenchmark%%.*}:%3R" # %3U %3S"
        echo -n "|-- Executing the script with bash..."
        { time /bin/bash "$script_to_execute" > $seq_output ; } \
            < "$stdin_redir" 2>> "${seq_time}"
        echo "   exited with $?"
        tail -n1 ${seq_time} >> ${results_time_bash}
        for conf in "${configurations[@]}"; do
            for n_in in "${n_inputs[@]}"; do
                echo "|-- Executing with pash --width ${n_in} ${conf}..."
                export pash_time="${test_results_dir}/${microbenchmark}_${n_in}_distr_$(echo ${conf} | tr -d ' ').time"
                export pash_output="${intermediary_dir}/${microbenchmark}_${n_in}_pash_output"
                export script_conf=${microbenchmark}_${n_in}
                echo '' > "${pash_time}"
                # do we need to write the PaSh output ?
                cat $stdin_redir |
                    execute_pash_and_check_diff -d $PASH_LOG $assert_correctness ${conf} --width "${n_in}" --output_time $script_to_execute                 
                tail -n1 "${pash_time}" >> "${results_time_pash}_${n_in}"
            done
        done
    done
}

execute_tests "" "${script_microbenchmarks[@]}"
execute_tests "--assert_all_regions_parallelizable" "${pipeline_microbenchmarks[@]}"

#cat ${results_time} | sed 's/,/./' > /tmp/a
#cat /tmp/a | sed 's/@/,/' > ${results_time}


if type lsb_release >/dev/null 2>&1 ; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# now do different things depending on distro
case "$distro" in
    freebsd*)  
        # change sed to gsed
        sed () {
            gsed $@
        }
        ;;
    *)
        ;;
esac

echo "group,Bash,Pash2,Pash8" > ${results_time}
paste -d'@' $test_results_dir/results.time_*  | sed 's\,\.\g' | sed 's\:\,\g' | sed 's\@\,\g' >> ${results_time}

#echo "Below follow the identical outputs:"
#grep "are identical" "$test_results_dir"/result_status | awk '{print $1}'

echo "Below follow the non-identical outputs:"     
grep "are not identical" "$test_results_dir"/result_status | awk '{print $1}'

TOTAL_TESTS=$(cat "$test_results_dir"/result_status | wc -l)
PASSED_TESTS=$(grep -c "are identical" "$test_results_dir"/result_status)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
