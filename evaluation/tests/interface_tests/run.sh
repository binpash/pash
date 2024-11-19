#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# time: print real in seconds, to simplify parsing

bash="bash"
pash="$PASH_TOP/pa.sh --profile_driven"

output_dir="$PASH_TOP/evaluation/tests/interface_tests/output"
rm -rf "$output_dir"
mkdir -p "$output_dir"

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
        echo "$test are not identical" >> $output_dir/result_status
        echo -e '\t\tFAIL'
        return 1
    else
        echo "$test are identical" >> $output_dir/result_status
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

## This test checks the behavior of the shell when -v is set.
##
## We do not care about what exactly is in the log, but rather about whether it contains something.
test12()
{
    local shell=$1
    $shell set-v.sh 2> set-v.log
    grep "echo hello" < set-v.log > /dev/null
    echo $?
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

test15()
{
    local shell=$1
    $shell readonly.sh 
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

test_set()
{
    local shell=$1
    $shell set.sh > /dev/null 2> set.sh.out
    cat set.sh.out
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

test_new_line_in_var()
{
    local shell=$1
    $shell newline_in_var.sh
}

test_cmd_sbst()
{
    local shell=$1
    `true; false; true; $shell cmd_sbst_subscript.sh`
}

test_cmd_sbst2()
{
    local shell=$1
    $shell cmd_sbst.sh
}

test_exec_redirections()
{
    local shell=$1
    $shell exec-redirections.sh
    cat exec-redirections.out
    [ -s exec-redirections.err ]
}

test_cat_hyphen()
{
    local shell=$1
    echo popo > /tmp/cat.in
    echo pipis >> /tmp/cat.in
    echo papa >> /tmp/cat.in
    echo "pipi" | $shell -c 'cat /tmp/cat.in - /tmp/cat.in | grep pipi'
}

test_trap()
{
    local shell=$1
    $shell trap.sh
}

test_set-dash()
{
    local shell=$1
    $shell -v -x set-dash-v-x.sh one two three four five > set-dash-v-x.out 2> set-dash-v-x.err
    grep 'echo hello' set-dash-v-x.err | wc -l
}

test_cat_redir_fail()
{
    local shell=$1
    $shell cat-redir-fail.sh
}

test_umask()
{
    local shell=$1
    umask u=rw,go-rwx
    $shell -c 'echo hi'
}

test_expand_u()
{
    local shell=$1
    $shell expand-u.sh
}

test_expand_u_positional()
{
    local shell=$1
    $shell expand-u-positional.sh
}

test_quoting()
{
    local shell=$1
    echo "ababa" | $shell -c 'tr -dc abc'
}

test_var_assgn_default()
{
    local shell=$1
    $shell var_assgn.sh
}

test_exclam()
{
    local shell=$1
    $shell test-exclam.sh
}

test_redir_var_test()
{
    local shell=$1
    $shell redir-var-test.sh
}

test_star()
{
    local shell=$1
    $shell test-star.sh foo '*' baz 'hi michael' "abc
     dfg"
}

test_env_vars()
{
    local shell=$1
    rm -f tmp1.txt tmp2.txt
    $shell env_vars.sh
    diff tmp1.txt tmp2.txt
}

test_redir_dup()
{
    local shell=$1
    $shell redir-dup.sh
}

test_IFS()
{
    local shell=$1
    $shell test-IFS.sh
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
    run_test test_set
    run_test test_set_e
    run_test test_redirect
    run_test test_unparsing
    run_test test_set_e_2
    run_test test_set_e_3
    run_test test_new_line_in_var
    run_test test_cmd_sbst
    run_test test_cmd_sbst2
    run_test test_exec_redirections
    run_test test_cat_hyphen
    run_test test_trap
    run_test test_set-dash
    run_test test_cat_redir_fail
    run_test test_umask
    run_test test_expand_u
    run_test test_expand_u_positional
    run_test test_quoting
    run_test test_var_assgn_default
    run_test test_exclam
    run_test test_redir_var_test
    run_test test_star
    run_test test_env_vars
    run_test test_redir_dup
    run_test test_IFS
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
grep "are identical" "$output_dir"/result_status | awk '{print $1}'

echo "Below follow the non-identical outputs:"     
grep "are not identical" "$output_dir"/result_status | awk '{print $1}'

TOTAL_TESTS=$(cat "$output_dir"/result_status | wc -l)
PASSED_TESTS=$(grep -c "are identical" "$output_dir"/result_status)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
