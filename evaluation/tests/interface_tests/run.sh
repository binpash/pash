#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# time: print real in seconds, to simplify parsing

bash="bash"
pash="$PASH_TOP/pa.sh"

output_dir="$PASH_TOP/evaluation/tests/interface_tests/output"
mkdir -p "$output_dir"

rm -f  $output_dir/results.time_bash
rm -f  $output_dir/results.time_pash

run_test()
{
    local test=$1
    local pash_specific_args=$2

    if [ "$(type -t $test)" != "function" ]; then
        echo "$test is not a function!   FAIL"
        return 1
    fi

    echo -n "Running $test..."
    TIMEFORMAT="${test%%.*}:%3R" # %3U %3S"
    { time $test "bash" > "$output_dir/$test.bash.out"; } 2>> $output_dir/results.time_bash
    test_bash_ec=$?
    TIMEFORMAT="%3R" # %3U %3S"
    { time $test "$pash" "$pash_specific_args" > "$output_dir/$test.pash.out"; } 2>> $output_dir/results.time_pash
    test_pash_ec=$?
    diff "$output_dir/$test.bash.out" "$output_dir/$test.pash.out"
    test_diff_ec=$?

    ## Check if the two exit codes are both success or both error
    { [ $test_bash_ec -eq 0 ] && [ $test_pash_ec -eq 0 ]; } || { [ $test_bash_ec -ne 0 ] && [ $test_pash_ec -ne 0 ]; }
    test_ec=$?
    if [ $test_diff_ec -ne 0 ]; then
        echo -n "$test output mismatch "
    fi
    if [ $test_ec -ne 0 ]; then
        echo -n "$test exit code mismatch "
    fi
    if [ $test_diff_ec -ne 0 ] || [ $test_ec -ne 0 ]; then
        echo "are not identical" > $output_dir/${test}_distr.time
        echo -e '\t\tFAIL'
        return 1
    else
        echo "are identical" > $output_dir/${test}_distr.time
        echo -e '\t\tOK'
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

test7()
{
    local shell=$1
    $shell args_with_spaces.sh "hello there" "hi friend"
}

test8()
{
    local shell=$1
    $shell -c 'shift 5; echo $#' pash 2 3 4
}

test9()
{
    local shell=$1
    $shell tilde.sh
}

test10()
{
    local shell=$1
    $shell exit_error.sh
}

test11()
{
    local shell=$1
    $shell incomplete-arith.sh
}

test12()
{
    local shell=$1
    $shell set-v.sh 2> set-v.log
    grep "echo hello" < set-v.log
}

test13()
{
    local shell=$1
    $shell -a -c 'var=1; env' > test13.log
    grep -e "^var=1" < test13.log
}

## Checks if +a is parsed correctly
test14()
{
    local shell=$1
    $shell +a readonly.sh
}

## Checks interactivity
##
## TODO: Make the interactivity script more elaborate (variable dependencies)
test15()
{
    local shell=$1
    $shell < readonly.sh 
}

test16()
{
    local shell=$1
    $shell heredoc2.sh
}

test17()
{
    local shell=$1
    $shell star-escape.sh
}

test18()
{
    local shell=$1
    $shell escape-madness.sh
}

test_set_e()
{
    local shell=$1
    $shell set-e.sh
}

test_set_e_2()
{
    local shell=$1
    $shell set-e-2.sh
}

test_set_e_3()
{
    local shell=$1
    $shell set-e-3.sh
}

test_redirect()
{
    local shell=$1
    $shell redirect_wrapper.sh "$shell"
}

test_unparsing()
{
    local shell=$1
    $shell unparsing-special-chars.sh
}

test_unset_exit_codes()
{
    local shell=$1
    local pash_specific_args=$2
    $shell $pash_specific_args ec-issue-127.sh
    echo $?
}

## We run all tests composed with && to exit on the first that fails
if [ "$#" -eq 0 ]; then
    run_test test1
    run_test test2
    run_test test3
    run_test test4
    run_test test5
    run_test test6
    # run_test test7
    run_test test8
    run_test test9
    run_test test10
    run_test test11
    run_test test12
    run_test test13
    run_test test14
    run_test test15
    run_test test16
    run_test test17
    run_test test18
    run_test test_set_e
    run_test test_redirect
    run_test test_unparsing
    run_test test_set_e_2
    run_test test_set_e_3
    run_test test_unset_exit_codes "--batch"
else
    for testname in $@
    do
        run_test "$testname"
    done
fi

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

echo "group,Bash,Pash-DRY_COMP" > $output_dir/results.time
paste $output_dir/results.time_*  | sed 's\,\.\g' | sed 's\:\,\g' | sed 's/\t/,/' >> $output_dir/results.time

echo "Below follow the identical outputs:"
grep --files-with-match "are identical" "$output_dir"/*_distr*.time

echo "Below follow the non-identical outputs:"
grep -L "are identical" "$output_dir"/*_distr*.time

TOTAL_TESTS=$(ls -la "$output_dir"/*_distr*.time | wc -l)
PASSED_TESTS=$(grep --files-with-match "are identical" "$output_dir"/*_distr*.time | wc -l)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
