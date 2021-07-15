#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# time: print real in seconds, to simplify parsing

## TODO: Make the compiler work too.
bash="bash"
pash="$PASH_TOP/pa.sh --dry_run_compiler"

output_dir="$PASH_TOP/evaluation/tests/interface_tests/output"
mkdir -p "$output_dir"

rm -f  $output_dir/results.time_bash
rm -f  $output_dir/results.time_pash
run_test()
{
    local test=$1

    if [ "$(type -t $test)" != "function" ]; then
        echo "$test is not a function!   FAIL"
        return 1
    fi

    echo -n "Running $test..."
    TIMEFORMAT="${test%%.*}:%3R" # %3U %3S"
    { time $test "bash" > "$output_dir/$test.bash.out"; } 2>> $output_dir/results.time_bash
    test_bash_ec=$?
    TIMEFORMAT="%3R" # %3U %3S"
    { time $test "$pash" > "$output_dir/$test.pash.out"; } 2>> $output_dir/results.time_pash
    test_pash_ec=$?
    diff "$output_dir/$test.bash.out" "$output_dir/$test.pash.out"
    test_diff_ec=$?

    if [ $test_diff_ec -ne 0 ]; then
        echo -n "$test output mismatch"
    fi
    if [ $test_bash_ec -ne $test_pash_ec ]; then
        echo -n "$test exit code mismatch"
    fi
    if [ $test_diff_ec -ne 0 ] || [ $test_bash_ec -ne $test_pash_ec ]; then
        echo "are not identical" > $output_dir/${test}_distr.time
        echo '   FAIL'
        return 1
    else
        echo "are identical" > $output_dir/${test}_distr.time
        echo '   OK'
        return 0
    fi
}

test1()
{
    local shell=$1
    $shell echo_args.sh 1 2 3
}

## Tests -c
test2()
{
    local shell=$1
    $shell -c 'echo $2 $1' pash 2 3
}

## Tests -c and the assignment to $0
test3()
{
    local shell=$1
    $shell -c 'echo $0 $2 $1' pash 2 3
}

test4()
{
    local shell=$1
    $shell -c 'shift; echo $1 $2' pash 2 3 4 5
}

test5()
{
    local shell=$1
    $shell heredoc1.sh
}

test6()
{
    local shell=$1
    $shell loop1.sh
}

## We run all tests composed with && to exit on the first that fails
if [ "$#" -eq 0 ]; then
    run_test test1 &&
    run_test test2 &&
    run_test test3 &&
    run_test test4 &&
    run_test test5 &&
    run_test test6

    # run_test test3  # This is commented out at the moment because it doesn't suceed
else
    for testname in $@
    do
        run_test "$testname"
    done
fi

echo "group,Bash,Pash-DRY_COMP" > $output_dir/results.time
paste $output_dir/results.time_*  | sed 's\,\.\g' | sed 's\:\,\g' | sed 's/\t/,/' >> $output_dir/results.time

echo "Below follow the identical outputs:"
grep --files-with-match "are identical" "$output_dir"/*_distr*.time

echo "Below follow the non-identical outputs:"
grep -L "are identical" "$output_dir"/*_distr*.time

TOTAL_TESTS=$(ls -la "$output_dir"/*_distr*.time | wc -l)
PASSED_TESTS=$(grep --files-with-match "are identical" "$output_dir"/*_distr*.time | wc -l)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
