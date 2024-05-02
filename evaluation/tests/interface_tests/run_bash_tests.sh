#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# time: print real in seconds, to simplify parsing

# are we running in bash mode?
test_mode=${1:-dash}
if [ "$test_mode" = "bash" ]; then
    echo "Running in bash mode"
fi

bash="bash"
if [ "$test_mode" = "bash" ]; then
    echo "Running in bash mode confirmed"
    pash="$PASH_TOP/pa.sh --parallel_pipelines --profile_driven --bash"
else
    pash="$PASH_TOP/pa.sh --parallel_pipelines --profile_driven"
fi

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

test_jobs1.sub()
{
    local shell=$1
    $shell bash_tests/jobs1.sub
}

test_exec8.sub()
{
    local shell=$1
    $shell bash_tests/exec8.sub
}

test_case1.sub()
{
    local shell=$1
    $shell bash_tests/case1.sub
}

test_nameref17.sub()
{
    local shell=$1
    $shell bash_tests/nameref17.sub
}

test_lastpipe1.sub()
{
    local shell=$1
    $shell bash_tests/lastpipe1.sub
}

test_dollar-at7.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at7.sub
}

test_trap4.sub()
{
    local shell=$1
    $shell bash_tests/trap4.sub
}

test_errors9.sub()
{
    local shell=$1
    $shell bash_tests/errors9.sub
}

test_dollar-star3.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star3.sub
}

test_assoc1.sub()
{
    local shell=$1
    $shell bash_tests/assoc1.sub
}

test_alias4.sub()
{
    local shell=$1
    $shell bash_tests/alias4.sub
}

test_varenv6.sub()
{
    local shell=$1
    $shell bash_tests/varenv6.sub
}

test_builtins2.sub()
{
    local shell=$1
    $shell bash_tests/builtins2.sub
}

test_array25.sub()
{
    local shell=$1
    $shell bash_tests/array25.sub
}

test_nquote1.sub()
{
    local shell=$1
    $shell bash_tests/nquote1.sub
}

test_array19.sub()
{
    local shell=$1
    $shell bash_tests/array19.sub
}

test_array18.sub()
{
    local shell=$1
    $shell bash_tests/array18.sub
}

test_redir1.sub()
{
    local shell=$1
    $shell bash_tests/redir1.sub
}

test_array24.sub()
{
    local shell=$1
    $shell bash_tests/array24.sub
}

test_array30.sub()
{
    local shell=$1
    $shell bash_tests/array30.sub
}

test_builtins3.sub()
{
    local shell=$1
    $shell bash_tests/builtins3.sub
}

test_varenv7.sub()
{
    local shell=$1
    $shell bash_tests/varenv7.sub
}

test_alias5.sub()
{
    local shell=$1
    $shell bash_tests/alias5.sub
}

test_dollar-star2.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star2.sub
}

test_nameref1.sub()
{
    local shell=$1
    $shell bash_tests/nameref1.sub
}

test_errors8.sub()
{
    local shell=$1
    $shell bash_tests/errors8.sub
}

test_trap5.sub()
{
    local shell=$1
    $shell bash_tests/trap5.sub
}

test_dollar-at6.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at6.sub
}

test_nameref16.sub()
{
    local shell=$1
    $shell bash_tests/nameref16.sub
}

test_exec9.sub()
{
    local shell=$1
    $shell bash_tests/exec9.sub
}

test_heredoc1.sub()
{
    local shell=$1
    $shell bash_tests/heredoc1.sub
}

test_read1.sub()
{
    local shell=$1
    $shell bash_tests/read1.sub
}

test_comsub-eof.tests()
{
    local shell=$1
    $shell bash_tests/comsub-eof.tests
}

test_jobs2.sub()
{
    local shell=$1
    $shell bash_tests/jobs2.sub
}

test_read3.sub()
{
    local shell=$1
    $shell bash_tests/read3.sub
}

test_heredoc3.sub()
{
    local shell=$1
    $shell bash_tests/heredoc3.sub
}

test_trap2a.sub()
{
    local shell=$1
    $shell bash_tests/trap2a.sub
}

test_case2.sub()
{
    local shell=$1
    $shell bash_tests/case2.sub
}

test_nameref14.sub()
{
    local shell=$1
    $shell bash_tests/nameref14.sub
}

test_lastpipe2.sub()
{
    local shell=$1
    $shell bash_tests/lastpipe2.sub
}

test_dollar-at4.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at4.sub
}

test_nameref3.sub()
{
    local shell=$1
    $shell bash_tests/nameref3.sub
}

test_quotearray.tests()
{
    local shell=$1
    $shell bash_tests/quotearray.tests
}

test_getopts8.sub()
{
    local shell=$1
    $shell bash_tests/getopts8.sub
}

test_assoc2.sub()
{
    local shell=$1
    $shell bash_tests/assoc2.sub
}

test_builtins1.sub()
{
    local shell=$1
    $shell bash_tests/builtins1.sub
}

test_dollar-at-star8.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star8.sub
}

test_varenv5.sub()
{
    local shell=$1
    $shell bash_tests/varenv5.sub
}

test_type4.sub()
{
    local shell=$1
    $shell bash_tests/type4.sub
}

test_array26.sub()
{
    local shell=$1
    $shell bash_tests/array26.sub
}

test_nquote2.sub()
{
    local shell=$1
    $shell bash_tests/nquote2.sub
}

test_redir3.sub()
{
    local shell=$1
    echo hi1 hi2 hi3 hi4 | $shell bash_tests/redir3.sub
}

test_redir2.sub()
{
    local shell=$1
    $shell bash_tests/redir2.sub
}

test_nquote3.sub()
{
    local shell=$1
    $shell bash_tests/nquote3.sub
}

test_array27.sub()
{
    local shell=$1
    $shell bash_tests/array27.sub
}

test_varenv4.sub()
{
    local shell=$1
    $shell bash_tests/varenv4.sub
}

test_dollar-at-star9.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star9.sub
}

test_assoc3.sub()
{
    local shell=$1
    $shell bash_tests/assoc3.sub
}

test_alias6.sub()
{
    local shell=$1
    $shell bash_tests/alias6.sub
}

test_posix2syntax.sub()
{
    local shell=$1
    $shell bash_tests/posix2syntax.sub
}

test_dollar-star1.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star1.sub
}

test_getopts9.sub()
{
    local shell=$1
    $shell bash_tests/getopts9.sub
}

test_dstack2.tests()
{
    local shell=$1
    $shell bash_tests/dstack2.tests
}

test_trap6.sub()
{
    local shell=$1
    $shell bash_tests/trap6.sub
}

test_nameref2.sub()
{
    local shell=$1
    $shell bash_tests/nameref2.sub
}

test_dollar-at5.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at5.sub
}

test_arith.tests()
{
    local shell=$1
    $shell bash_tests/arith.tests
}

test_lastpipe3.sub()
{
    local shell=$1
    $shell bash_tests/lastpipe3.sub
}

test_nameref15.sub()
{
    local shell=$1
    $shell bash_tests/nameref15.sub
}

test_case3.sub()
{
    local shell=$1
    $shell bash_tests/case3.sub
}

test_procsub.tests()
{
    local shell=$1
    $shell bash_tests/procsub.tests
}

test_heredoc2.sub()
{
    local shell=$1
    $shell bash_tests/heredoc2.sub
}

test_read2.sub()
{
    local shell=$1
    $shell bash_tests/read2.sub
}

test_jobs3.sub()
{
    local shell=$1
    $shell bash_tests/jobs3.sub
}

test_heredoc6.sub()
{
    local shell=$1
    $shell bash_tests/heredoc6.sub
}

test_jobs7.sub()
{
    local shell=$1
    $shell bash_tests/jobs7.sub
}

test_read6.sub()
{
    local shell=$1
    $shell bash_tests/read6.sub
}

test_shopt.tests()
{
    local shell=$1
    $shell bash_tests/shopt.tests
}

test_comsub.tests()
{
    local shell=$1
    $shell bash_tests/comsub.tests
}

test_dollar-star10.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star10.sub
}

test_dollar-at1.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at1.sub
}

test_appendop1.sub()
{
    local shell=$1
    $shell bash_tests/appendop1.sub
}

test_nameref11.sub()
{
    local shell=$1
    $shell bash_tests/nameref11.sub
}

test_strip.tests()
{
    local shell=$1
    $shell bash_tests/strip.tests
}

test_nameref6.sub()
{
    local shell=$1
    $shell bash_tests/nameref6.sub
}

test_trap2.sub()
{
    local shell=$1
    $shell bash_tests/trap2.sub
}

test_alias2.sub()
{
    local shell=$1
    $shell bash_tests/alias2.sub
}

test_assoc7.sub()
{
    local shell=$1
    $shell bash_tests/assoc7.sub
}

test_dollar-star5.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star5.sub
}

test_builtins4.sub()
{
    local shell=$1
    $shell bash_tests/builtins4.sub
}

test_new-exp9.sub()
{
    local shell=$1
    $shell bash_tests/new-exp9.sub
}

test_redir6.sub()
{
    local shell=$1
    $shell bash_tests/redir6.sub
}

test_array23.sub()
{
    local shell=$1
    $shell bash_tests/array23.sub
}

test_type1.sub()
{
    local shell=$1
    $shell bash_tests/type1.sub
}

test_array22.sub()
{
    local shell=$1
    $shell bash_tests/array22.sub
}

test_redir7.sub()
{
    local shell=$1
    $shell bash_tests/redir7.sub
}

test_intl1.sub()
{
    local shell=$1
    $shell bash_tests/intl1.sub
}

test_varenv1.sub()
{
    local shell=$1
    $shell bash_tests/varenv1.sub
}

test_builtins5.sub()
{
    local shell=$1
    $shell bash_tests/builtins5.sub
}

test_new-exp8.sub()
{
    local shell=$1
    $shell bash_tests/new-exp8.sub
}

test_dollar-star4.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star4.sub
}

test_assoc6.sub()
{
    local shell=$1
    $shell bash_tests/assoc6.sub
}

test_alias3.sub()
{
    local shell=$1
    $shell bash_tests/alias3.sub
}

test_set-x1.sub()
{
    local shell=$1
    $shell bash_tests/set-x1.sub
}

test_trap3.sub()
{
    local shell=$1
    $shell bash_tests/trap3.sub
}

test_nameref7.sub()
{
    local shell=$1
    $shell bash_tests/nameref7.sub
}

test_history.tests()
{
    local shell=$1
    $shell bash_tests/history.tests
}

test_nameref10.sub()
{
    local shell=$1
    $shell bash_tests/nameref10.sub
}

test_parser1.sub()
{
    local shell=$1
    $shell bash_tests/parser1.sub
}

test_jobs.tests()
{
    local shell=$1
    $shell bash_tests/jobs.tests
}

test_posixexp8.sub()
{
    local shell=$1
    $shell bash_tests/posixexp8.sub
}

test_read7.sub()
{
    local shell=$1
    $shell bash_tests/read7.sub
}

test_jobs6.sub()
{
    local shell=$1
    $shell bash_tests/jobs6.sub
}

test_heredoc7.sub()
{
    local shell=$1
    $shell bash_tests/heredoc7.sub
}

test_heredoc5.sub()
{
    local shell=$1
    $shell bash_tests/heredoc5.sub
}

test_read5.sub()
{
    local shell=$1
    $shell bash_tests/read5.sub
}

test_jobs4.sub()
{
    local shell=$1
    $shell bash_tests/jobs4.sub
}

test_case4.sub()
{
    local shell=$1
    $shell bash_tests/case4.sub
}

test_dollar-at2.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at2.sub
}

test_nameref.tests()
{
    local shell=$1
    $shell bash_tests/nameref.tests
}

test_nameref12.sub()
{
    local shell=$1
    $shell bash_tests/nameref12.sub
}

test_appendop2.sub()
{
    local shell=$1
    $shell bash_tests/appendop2.sub
}

test_trap1.sub()
{
    local shell=$1
    $shell bash_tests/trap1.sub
}

test_nameref5.sub()
{
    local shell=$1
    $shell bash_tests/nameref5.sub
}

test_assoc4.sub()
{
    local shell=$1
    $shell bash_tests/assoc4.sub
}

test_tilde2.tests()
{
    local shell=$1
    $shell bash_tests/tilde2.tests
}

test_vredir.tests()
{
    local shell=$1
    $shell bash_tests/vredir.tests
}

test_alias1.sub()
{
    local shell=$1
    $shell bash_tests/alias1.sub
}

test_mapfile.tests()
{
    local shell=$1
    $shell bash_tests/mapfile.tests
}

test_dollar-star6.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star6.sub
}

test_parser.tests()
{
    local shell=$1
    $shell bash_tests/parser.tests
}

test_varenv3.sub()
{
    local shell=$1
    $shell bash_tests/varenv3.sub
}

test_builtins7.sub()
{
    local shell=$1
    $shell bash_tests/builtins7.sub
}

test_intl3.sub()
{
    local shell=$1
    $shell bash_tests/intl3.sub
}

test_redir5.sub()
{
    local shell=$1
    $shell bash_tests/redir5.sub
}

test_errors.tests()
{
    local shell=$1
    $shell bash_tests/errors.tests
}

test_type2.sub()
{
    local shell=$1
    $shell bash_tests/type2.sub
}

test_array20.sub()
{
    local shell=$1
    $shell bash_tests/array20.sub
}

test_nquote4.sub()
{
    local shell=$1
    $shell bash_tests/nquote4.sub
}

test_type3.sub()
{
    local shell=$1
    $shell bash_tests/type3.sub
}

test_nquote5.sub()
{
    local shell=$1
    $shell bash_tests/nquote5.sub
}

test_array21.sub()
{
    local shell=$1
    $shell bash_tests/array21.sub
}

test_redir4.sub()
{
    local shell=$1
    $shell bash_tests/redir4.sub
}

test_intl2.sub()
{
    local shell=$1
    $shell bash_tests/intl2.sub
}

test_herestr1.sub()
{
    local shell=$1
    $shell bash_tests/herestr1.sub
}

test_builtins6.sub()
{
    local shell=$1
    $shell bash_tests/builtins6.sub
}

test_varenv2.sub()
{
    local shell=$1
    $shell bash_tests/varenv2.sub
}

test_dollar-star7.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star7.sub
}

test_assoc5.sub()
{
    local shell=$1
    $shell bash_tests/assoc5.sub
}

test_nameref4.sub()
{
    local shell=$1
    $shell bash_tests/nameref4.sub
}

test_nameref13.sub()
{
    local shell=$1
    $shell bash_tests/nameref13.sub
}

test_dollar-at3.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at3.sub
}

test_extglob3.tests()
{
    local shell=$1
    $shell bash_tests/extglob3.tests
}

test_jobs5.sub()
{
    local shell=$1
    $shell bash_tests/jobs5.sub
}

test_read4.sub()
{
    local shell=$1
    $shell bash_tests/read4.sub
}

test_heredoc4.sub()
{
    local shell=$1
    $shell bash_tests/heredoc4.sub
}

test_dbg-support3.sub()
{
    local shell=$1
    $shell bash_tests/dbg-support3.sub
}

test_globstar2.sub()
{
    local shell=$1
    $shell bash_tests/globstar2.sub
}

test_arith4.sub()
{
    local shell=$1
    $shell bash_tests/arith4.sub
}

test_vredir4.sub()
{
    local shell=$1
    $shell bash_tests/vredir4.sub
}

test_exp6.sub()
{
    local shell=$1
    $shell bash_tests/exp6.sub
}

test_comsub-eof3.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof3.sub
}

test_varenv.tests()
{
    local shell=$1
    $shell bash_tests/varenv.tests
}

test_rhs-exp.tests()
{
    local shell=$1
    $shell bash_tests/rhs-exp.tests
}

test_varenv12.sub()
{
    local shell=$1
    $shell bash_tests/varenv12.sub
}

test_attr1.sub()
{
    local shell=$1
    $shell bash_tests/attr1.sub
}

test_glob1.sub()
{
    local shell=$1
    $shell bash_tests/glob1.sub
}

test_glob.tests()
{
    local shell=$1
    $shell bash_tests/glob.tests
}

test_set-e3a.sub()
{
    local shell=$1
    $shell bash_tests/set-e3a.sub
}

test_varenv13.sub()
{
    local shell=$1
    $shell bash_tests/varenv13.sub
}

test_comsub-eof2.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof2.sub
}

test_exp7.sub()
{
    local shell=$1
    $shell bash_tests/exp7.sub
}

test_braces.tests()
{
    local shell=$1
    $shell bash_tests/braces.tests
}

test_vredir5.sub()
{
    local shell=$1
    $shell bash_tests/vredir5.sub
}

test_globstar3.sub()
{
    local shell=$1
    $shell bash_tests/globstar3.sub
}

test_arith5.sub()
{
    local shell=$1
    $shell bash_tests/arith5.sub
}

test_vredir7.sub()
{
    local shell=$1
    $shell bash_tests/vredir7.sub
}

test_arith7.sub()
{
    local shell=$1
    $shell bash_tests/arith7.sub
}

test_globstar1.sub()
{
    local shell=$1
    $shell bash_tests/globstar1.sub
}

test_comsub-eof0.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof0.sub
}

test_exp5.sub()
{
    local shell=$1
    $shell bash_tests/exp5.sub
}

test_varenv11.sub()
{
    local shell=$1
    $shell bash_tests/varenv11.sub
}

test_glob3.sub()
{
    local shell=$1
    $shell bash_tests/glob3.sub
}

test_glob10.sub()
{
    local shell=$1
    $shell bash_tests/glob10.sub
}

test_dbg-support2.tests()
{
    local shell=$1
    $shell bash_tests/dbg-support2.tests
}

test_quote4.sub()
{
    local shell=$1
    $shell bash_tests/quote4.sub
}

test_dynvar.tests()
{
    local shell=$1
    $shell bash_tests/dynvar.tests
}

test_case.tests()
{
    local shell=$1
    $shell bash_tests/case.tests
}

test_attr2.sub()
{
    local shell=$1
    $shell bash_tests/attr2.sub
}

test_extglob1a.sub()
{
    local shell=$1
    $shell bash_tests/extglob1a.sub
}

test_glob2.sub()
{
    local shell=$1
    $shell bash_tests/glob2.sub
}

test_varenv10.sub()
{
    local shell=$1
    $shell bash_tests/varenv10.sub
}

test_exp4.sub()
{
    local shell=$1
    $shell bash_tests/exp4.sub
}

test_comsub-eof1.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof1.sub
}

test_arith6.sub()
{
    local shell=$1
    $shell bash_tests/arith6.sub
}

test_vredir6.sub()
{
    local shell=$1
    $shell bash_tests/vredir6.sub
}

test_assoc18.sub()
{
    local shell=$1
    $shell bash_tests/assoc18.sub
}

test_shopt1.sub()
{
    local shell=$1
    $shell bash_tests/shopt1.sub
}

test_array8.sub()
{
    local shell=$1
    $shell bash_tests/array8.sub
}

test_vredir2.sub()
{
    local shell=$1
    $shell bash_tests/vredir2.sub
}

test_arith2.sub()
{
    local shell=$1
    $shell bash_tests/arith2.sub
}

test_varenv14.sub()
{
    local shell=$1
    $shell bash_tests/varenv14.sub
}

test_comsub-eof5.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof5.sub
}

test_rsh1.sub()
{
    local shell=$1
    $shell bash_tests/rsh1.sub
}

test_glob6.sub()
{
    local shell=$1
    $shell bash_tests/glob6.sub
}

test_arith-for.tests()
{
    local shell=$1
    $shell bash_tests/arith-for.tests
}

test_quote1.sub()
{
    local shell=$1
    $shell bash_tests/quote1.sub
}

test_read.tests()
{
    local shell=$1
    $shell bash_tests/read.tests
}

test_cond-regexp1.sub()
{
    local shell=$1
    $shell bash_tests/cond-regexp1.sub
}

test_extglob2.tests()
{
    local shell=$1
    $shell bash_tests/extglob2.tests
}

test_appendop.tests()
{
    local shell=$1
    $shell bash_tests/appendop.tests
}

test_iquote1.sub()
{
    local shell=$1
    $shell bash_tests/iquote1.sub
}

test_exportfunc1.sub()
{
    local shell=$1
    $shell bash_tests/exportfunc1.sub
}

test_extglob.tests()
{
    local shell=$1
    $shell bash_tests/extglob.tests
}

test_glob7.sub()
{
    local shell=$1
    $shell bash_tests/glob7.sub
}

test_dbg-support.tests()
{
    local shell=$1
    $shell bash_tests/dbg-support.tests
}

test_exp1.sub()
{
    local shell=$1
    $shell bash_tests/exp1.sub
}

test_comsub-eof4.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof4.sub
}

test_rsh.tests()
{
    local shell=$1
    $shell bash_tests/rsh.tests
}

test_varenv15.sub()
{
    local shell=$1
    $shell bash_tests/varenv15.sub
}

test_func.tests()
{
    local shell=$1
    $shell bash_tests/func.tests
}

test_arith3.sub()
{
    local shell=$1
    $shell bash_tests/arith3.sub
}

test_vredir3.sub()
{
    local shell=$1
    $shell bash_tests/vredir3.sub
}

test_rhs-exp1.sub()
{
    local shell=$1
    $shell bash_tests/rhs-exp1.sub
}

test_array9.sub()
{
    local shell=$1
    $shell bash_tests/array9.sub
}

test_posixexp.tests()
{
    local shell=$1
    $shell bash_tests/posixexp.tests
}

test_arith1.sub()
{
    local shell=$1
    $shell bash_tests/arith1.sub
}

test_vredir1.sub()
{
    local shell=$1
    $shell bash_tests/vredir1.sub
}

test_cond.tests()
{
    local shell=$1
    $shell bash_tests/cond.tests
}

test_ifs.tests()
{
    local shell=$1
    $shell bash_tests/ifs.tests
}

test_varenv17.sub()
{
    local shell=$1
    $shell bash_tests/varenv17.sub
}

test_exp3.sub()
{
    local shell=$1
    $shell bash_tests/exp3.sub
}

test_comsub-eof6.sub()
{
    local shell=$1
    $shell bash_tests/comsub-eof6.sub
}

test_rsh2.sub()
{
    local shell=$1
    $shell bash_tests/rsh2.sub
}

test_glob5.sub()
{
    local shell=$1
    $shell bash_tests/glob5.sub
}

test_quote2.sub()
{
    local shell=$1
    $shell bash_tests/quote2.sub
}

test_lastpipe.tests()
{
    local shell=$1
    $shell bash_tests/lastpipe.tests
}

test_exportfunc3.sub()
{
    local shell=$1
    $shell bash_tests/exportfunc3.sub
}

test_cond-regexp2.sub()
{
    local shell=$1
    $shell bash_tests/cond-regexp2.sub
}

test_more-exp.tests()
{
    local shell=$1
    $shell bash_tests/more-exp.tests
}

test_cond-regexp3.sub()
{
    local shell=$1
    $shell bash_tests/cond-regexp3.sub
}

test_exportfunc2.sub()
{
    local shell=$1
    $shell bash_tests/exportfunc2.sub
}

test_quote3.sub()
{
    local shell=$1
    $shell bash_tests/quote3.sub
}

test_glob4.sub()
{
    local shell=$1
    $shell bash_tests/glob4.sub
}

test_exp2.sub()
{
    local shell=$1
    $shell bash_tests/exp2.sub
}

test_varenv16.sub()
{
    local shell=$1
    $shell bash_tests/varenv16.sub
}

test_assoc13.sub()
{
    local shell=$1
    $shell bash_tests/assoc13.sub
}

test_array7.sub()
{
    local shell=$1
    $shell bash_tests/array7.sub
}

test_exp.tests()
{
    local shell=$1
    $shell bash_tests/exp.tests
}

test_new-exp12.sub()
{
    local shell=$1
    $shell bash_tests/new-exp12.sub
}

test_new-exp.tests()
{
    local shell=$1
    $shell bash_tests/new-exp.tests
}

test_exp12.sub()
{
    local shell=$1
    $shell bash_tests/exp12.sub
}

test_comsub6.sub()
{
    local shell=$1
    $shell bash_tests/comsub6.sub
}

test_histexp6.sub()
{
    local shell=$1
    $shell bash_tests/histexp6.sub
}

test_getopts.tests()
{
    local shell=$1
    $shell bash_tests/getopts.tests
}

test_glob9.sub()
{
    local shell=$1
    $shell bash_tests/glob9.sub
}

test_input-line.sub()
{
    local shell=$1
    echo i | $shell bash_tests/input-line.sub
}

test_history3.sub()
{
    local shell=$1
    $shell bash_tests/history3.sub
}

test_func1.sub()
{
    local shell=$1
    $shell bash_tests/func1.sub
}

test_extglob4.sub()
{
    local shell=$1
    $shell bash_tests/extglob4.sub
}

test_dstack.tests()
{
    local shell=$1
    $shell bash_tests/dstack.tests
}

test_extglob5.sub()
{
    local shell=$1
    $shell bash_tests/extglob5.sub
}

test_history2.sub()
{
    local shell=$1
    $shell bash_tests/history2.sub
}

test_unicode1.sub()
{
    local shell=$1
    $shell bash_tests/unicode1.sub
}

test_glob8.sub()
{
    local shell=$1
    $shell bash_tests/glob8.sub
}

test_histexp7.sub()
{
    local shell=$1
    $shell bash_tests/histexp7.sub
}

test_exp13.sub()
{
    local shell=$1
    $shell bash_tests/exp13.sub
}

test_new-exp13.sub()
{
    local shell=$1
    $shell bash_tests/new-exp13.sub
}

test_assoc.tests()
{
    local shell=$1
    $shell bash_tests/assoc.tests
}

test_dbg-support.sub()
{
    local shell=$1
    $shell bash_tests/dbg-support.sub
}

test_array6.sub()
{
    local shell=$1
    $shell bash_tests/array6.sub
}

test_assoc12.sub()
{
    local shell=$1
    $shell bash_tests/assoc12.sub
}

test_assoc10.sub()
{
    local shell=$1
    $shell bash_tests/assoc10.sub
}

test_array4.sub()
{
    local shell=$1
    $shell bash_tests/array4.sub
}

test_posixpipe.tests()
{
    local shell=$1
    $shell bash_tests/posixpipe.tests
}

test_new-exp11.sub()
{
    local shell=$1
    $shell bash_tests/new-exp11.sub
}

test_exp11.sub()
{
    local shell=$1
    $shell bash_tests/exp11.sub
}

test_comsub5.sub()
{
    local shell=$1
    $shell bash_tests/comsub5.sub
}

test_varenv18.sub()
{
    local shell=$1
    $shell bash_tests/varenv18.sub
}

test_posix2.tests()
{
    local shell=$1
    $shell bash_tests/posix2.tests
}

test_histexp5.sub()
{
    local shell=$1
    $shell bash_tests/histexp5.sub
}

test_printf.tests()
{
    local shell=$1
    $shell bash_tests/printf.tests
}

test_unicode3.sub()
{
    local shell=$1
    $shell bash_tests/unicode3.sub
}

test_set-x.tests()
{
    local shell=$1
    $shell bash_tests/set-x.tests
}

test_extglob7.sub()
{
    local shell=$1
    $shell bash_tests/extglob7.sub
}

test_func2.sub()
{
    local shell=$1
    $shell bash_tests/func2.sub
}

test_comsub-posix.tests()
{
    local shell=$1
    $shell bash_tests/comsub-posix.tests
}

test_quote.tests()
{
    local shell=$1
    $shell bash_tests/quote.tests
}

test_extglob6.sub()
{
    local shell=$1
    $shell bash_tests/extglob6.sub
}

test_func3.sub()
{
    local shell=$1
    $shell bash_tests/func3.sub
}

test_history1.sub()
{
    local shell=$1
    $shell bash_tests/history1.sub
}

test_unicode2.sub()
{
    local shell=$1
    $shell bash_tests/unicode2.sub
}

test_tilde3.sub()
{
    local shell=$1
    $shell bash_tests/tilde3.sub
}

test_casemod.tests()
{
    local shell=$1
    $shell bash_tests/casemod.tests
}

test_histexp4.sub()
{
    local shell=$1
    $shell bash_tests/histexp4.sub
}

test_comsub4.sub()
{
    local shell=$1
    $shell bash_tests/comsub4.sub
}

test_getopts10.sub()
{
    local shell=$1
    $shell bash_tests/getopts10.sub
}

test_varenv19.sub()
{
    local shell=$1
    $shell bash_tests/varenv19.sub
}

test_exp10.sub()
{
    local shell=$1
    $shell bash_tests/exp10.sub
}

test_nquote5.tests()
{
    local shell=$1
    $shell bash_tests/nquote5.tests
}

test_new-exp10.sub()
{
    local shell=$1
    $shell bash_tests/new-exp10.sub
}

test_array5.sub()
{
    local shell=$1
    $shell bash_tests/array5.sub
}

test_test1.sub()
{
    local shell=$1
    $shell bash_tests/test1.sub
}

test_assoc11.sub()
{
    local shell=$1
    $shell bash_tests/assoc11.sub
}

test_array1.sub()
{
    local shell=$1
    $shell bash_tests/array1.sub
}

test_assoc15.sub()
{
    local shell=$1
    $shell bash_tests/assoc15.sub
}

test_new-exp14.sub()
{
    local shell=$1
    $shell bash_tests/new-exp14.sub
}

test_varenv21.sub()
{
    local shell=$1
    $shell bash_tests/varenv21.sub
}

test_exp9.sub()
{
    local shell=$1
    $shell bash_tests/exp9.sub
}

test_procsub2.sub()
{
    local shell=$1
    $shell bash_tests/procsub2.sub
}

test_histexp.tests()
{
    local shell=$1
    $shell bash_tests/histexp.tests
}

test_history5.sub()
{
    local shell=$1
    $shell bash_tests/history5.sub
}

test_nquote.tests()
{
    local shell=$1
    $shell bash_tests/nquote.tests
}

test_extglob2.sub()
{
    local shell=$1
    $shell bash_tests/extglob2.sub
}

test_nquote1.tests()
{
    local shell=$1
    $shell bash_tests/nquote1.tests
}

test_set-e1.sub()
{
    local shell=$1
    $shell bash_tests/set-e1.sub
}

test_extglob3.sub()
{
    local shell=$1
    $shell bash_tests/extglob3.sub
}

test_globstar.tests()
{
    local shell=$1
    $shell bash_tests/globstar.tests
}

test_history4.sub()
{
    local shell=$1
    $shell bash_tests/history4.sub
}

test_intl.tests()
{
    local shell=$1
    $shell bash_tests/intl.tests
}

test_comsub1.sub()
{
    local shell=$1
    $shell bash_tests/comsub1.sub
}

test_histexp1.sub()
{
    local shell=$1
    $shell bash_tests/histexp1.sub
}

test_varenv20.sub()
{
    local shell=$1
    $shell bash_tests/varenv20.sub
}

test_exp8.sub()
{
    local shell=$1
    $shell bash_tests/exp8.sub
}

test_new-exp15.sub()
{
    local shell=$1
    $shell bash_tests/new-exp15.sub
}

test_complete.tests()
{
    local shell=$1
    $shell bash_tests/complete.tests
}

test_assoc14.sub()
{
    local shell=$1
    $shell bash_tests/assoc14.sub
}

test_array2.sub()
{
    local shell=$1
    $shell bash_tests/array2.sub
}

test_assoc16.sub()
{
    local shell=$1
    $shell bash_tests/assoc16.sub
}

test_vredir8.sub()
{
    local shell=$1
    $shell bash_tests/vredir8.sub
}

test_cond-regexp.sub()
{
    local shell=$1
    $shell bash_tests/cond-regexp.sub
}

test_arith8.sub()
{
    local shell=$1
    $shell bash_tests/arith8.sub
}

test_set-e.tests()
{
    local shell=$1
    $shell bash_tests/set-e.tests
}

test_histexp3.sub()
{
    local shell=$1
    $shell bash_tests/histexp3.sub
}

test_varenv22.sub()
{
    local shell=$1
    $shell bash_tests/varenv22.sub
}

test_procsub1.sub()
{
    local shell=$1
    $shell bash_tests/procsub1.sub
}

test_posixexp2.tests()
{
    local shell=$1
    $shell bash_tests/posixexp2.tests
}

test_comsub3.sub()
{
    local shell=$1
    $shell bash_tests/comsub3.sub
}

test_history6.sub()
{
    local shell=$1
    $shell bash_tests/history6.sub
}

test_builtins.tests()
{
    local shell=$1
    $shell bash_tests/builtins.tests
}

test_extglob1.sub()
{
    local shell=$1
    $shell bash_tests/extglob1.sub
}

test_func4.sub()
{
    local shell=$1
    $shell bash_tests/func4.sub
}

test_alias.tests()
{
    local shell=$1
    $shell bash_tests/alias.tests
}

test_set-e3.sub()
{
    local shell=$1
    $shell bash_tests/set-e3.sub
}

test_set-e2.sub()
{
    local shell=$1
    $shell bash_tests/set-e2.sub
}

test_nquote3.tests()
{
    local shell=$1
    $shell bash_tests/nquote3.tests
}

test_comsub2.sub()
{
    local shell=$1
    $shell bash_tests/comsub2.sub
}

test_histexp2.sub()
{
    local shell=$1
    $shell bash_tests/histexp2.sub
}

test_new-exp16.sub()
{
    local shell=$1
    $shell bash_tests/new-exp16.sub
}

test_assoc17.sub()
{
    local shell=$1
    $shell bash_tests/assoc17.sub
}

test_array3.sub()
{
    local shell=$1
    $shell bash_tests/array3.sub
}

test_type.tests()
{
    local shell=$1
    $shell bash_tests/type.tests
}

test_exec1.sub()
{
    local shell=$1
    $shell bash_tests/exec1.sub
}

test_posixexp6.sub()
{
    local shell=$1
    $shell bash_tests/posixexp6.sub
}

test_redir11.sub()
{
    local shell=$1
    $shell bash_tests/redir11.sub
}

test_nameref22.sub()
{
    local shell=$1
    $shell bash_tests/nameref22.sub
}

test_getopts2.sub()
{
    local shell=$1
    $shell bash_tests/getopts2.sub
}

test_nameref9.sub()
{
    local shell=$1
    $shell bash_tests/nameref9.sub
}

test_assoc8.sub()
{
    local shell=$1
    $shell bash_tests/assoc8.sub
}

test_new-exp6.sub()
{
    local shell=$1
    $shell bash_tests/new-exp6.sub
}

test_dollar-at-star2.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star2.sub
}

test_exec12.sub()
{
    local shell=$1
    $shell bash_tests/exec12.sub
}

test_source1.sub()
{
    local shell=$1
    $shell bash_tests/source1.sub
}

test_array10.sub()
{
    local shell=$1
    $shell bash_tests/array10.sub
}

test_redir9.sub()
{
    local shell=$1
    $shell bash_tests/redir9.sub
}

test_quotearray2.sub()
{
    local shell=$1
    $shell bash_tests/quotearray2.sub
}

test_redir8.sub()
{
    local shell=$1
    $shell bash_tests/redir8.sub
}

test_quotearray3.sub()
{
    local shell=$1
    $shell bash_tests/quotearray3.sub
}

test_array11.sub()
{
    local shell=$1
    $shell bash_tests/array11.sub
}

test_dollar-at-star3.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star3.sub
}

test_exec13.sub()
{
    local shell=$1
    $shell bash_tests/exec13.sub
}

test_new-exp7.sub()
{
    local shell=$1
    $shell bash_tests/new-exp7.sub
}

test_tilde.tests()
{
    local shell=$1
    $shell bash_tests/tilde.tests
}

test_assoc9.sub()
{
    local shell=$1
    $shell bash_tests/assoc9.sub
}

test_nquote4.tests()
{
    local shell=$1
    $shell bash_tests/nquote4.tests
}

test_errors1.sub()
{
    local shell=$1
    $shell bash_tests/errors1.sub
}

test_nameref8.sub()
{
    local shell=$1
    $shell bash_tests/nameref8.sub
}

test_getopts3.sub()
{
    local shell=$1
    $shell bash_tests/getopts3.sub
}

test_ifs-posix.tests()
{
    local shell=$1
    $shell bash_tests/ifs-posix.tests
}

test_comsub-posix1.sub()
{
    local shell=$1
    $shell bash_tests/comsub-posix1.sub
}

test_nameref23.sub()
{
    local shell=$1
    $shell bash_tests/nameref23.sub
}

test_redir10.sub()
{
    local shell=$1
    $shell bash_tests/redir10.sub
}

test_printf4.sub()
{
    local shell=$1
    $shell bash_tests/printf4.sub
}

test_posixexp7.sub()
{
    local shell=$1
    $shell bash_tests/posixexp7.sub
}

test_read8.sub()
{
    local shell=$1
    $shell bash_tests/read8.sub
}

test_exec2.sub()
{
    local shell=$1
    $shell bash_tests/exec2.sub
}

test_posixexp5.sub()
{
    local shell=$1
    $shell bash_tests/posixexp5.sub
}

test_trap.tests()
{
    local shell=$1
    $shell bash_tests/trap.tests
}

test_redir.tests()
{
    local shell=$1
    $shell bash_tests/redir.tests
}

test_redir12.sub()
{
    local shell=$1
    $shell bash_tests/redir12.sub
}

test_nameref21.sub()
{
    local shell=$1
    $shell bash_tests/nameref21.sub
}

test_comsub-posix3.sub()
{
    local shell=$1
    $shell bash_tests/comsub-posix3.sub
}

test_getopts1.sub()
{
    local shell=$1
    $shell bash_tests/getopts1.sub
}

test_dollar-at-star10.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star10.sub
}

test_errors3.sub()
{
    local shell=$1
    $shell bash_tests/errors3.sub
}

test_dollar-star9.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star9.sub
}

test_exec11.sub()
{
    local shell=$1
    $shell bash_tests/exec11.sub
}

test_dollar-at-star1.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star1.sub
}

test_new-exp5.sub()
{
    local shell=$1
    $shell bash_tests/new-exp5.sub
}

test_source2.sub()
{
    local shell=$1
    $shell bash_tests/source2.sub
}

test_quotearray1.sub()
{
    local shell=$1
    $shell bash_tests/quotearray1.sub
}

test_array13.sub()
{
    local shell=$1
    $shell bash_tests/array13.sub
}

test_posixpat.tests()
{
    local shell=$1
    $shell bash_tests/posixpat.tests
}

test_array.tests()
{
    local shell=$1
    $shell bash_tests/array.tests
}

test_array12.sub()
{
    local shell=$1
    $shell bash_tests/array12.sub
}

test_iquote.tests()
{
    local shell=$1
    $shell bash_tests/iquote.tests
}

test_source3.sub()
{
    local shell=$1
    $shell bash_tests/source3.sub
}

test_new-exp4.sub()
{
    local shell=$1
    $shell bash_tests/new-exp4.sub
}

test_cprint.tests()
{
    local shell=$1
    $shell bash_tests/cprint.tests
}

test_exec10.sub()
{
    local shell=$1
    $shell bash_tests/exec10.sub
}

test_test.tests()
{
    local shell=$1
    $shell bash_tests/test.tests
}

test_dollar-star8.sub()
{
    local shell=$1
    $shell bash_tests/dollar-star8.sub
}

test_errors2.sub()
{
    local shell=$1
    $shell bash_tests/errors2.sub
}

test_dollar-at-star11.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star11.sub
}

test_comsub-posix2.sub()
{
    local shell=$1
    $shell bash_tests/comsub-posix2.sub
}

test_attr.tests()
{
    local shell=$1
    $shell bash_tests/attr.tests
}

test_nameref20.sub()
{
    local shell=$1
    $shell bash_tests/nameref20.sub
}

test_invert.tests()
{
    local shell=$1
    $shell bash_tests/invert.tests
}

test_posixexp4.sub()
{
    local shell=$1
    $shell bash_tests/posixexp4.sub
}

test_exec3.sub()
{
    local shell=$1
    $shell bash_tests/exec3.sub
}

test_exec7.sub()
{
    local shell=$1
    $shell bash_tests/exec7.sub
}

test_precedence.tests()
{
    local shell=$1
    $shell bash_tests/precedence.tests
}

test_printf3.sub()
{
    local shell=$1
    $shell bash_tests/printf3.sub
}

test_nameref18.sub()
{
    local shell=$1
    $shell bash_tests/nameref18.sub
}

test_comsub-posix6.sub()
{
    local shell=$1
    $shell bash_tests/comsub-posix6.sub
}

test_errors6.sub()
{
    local shell=$1
    $shell bash_tests/errors6.sub
}

test_getopts4.sub()
{
    local shell=$1
    $shell bash_tests/getopts4.sub
}

test_herestr.tests()
{
    local shell=$1
    $shell bash_tests/herestr.tests
}

test_source7.sub()
{
    local shell=$1
    $shell bash_tests/source7.sub
}

test_varenv9.sub()
{
    local shell=$1
    $shell bash_tests/varenv9.sub
}

test_dollar-at-star4.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star4.sub
}

test_exec14.sub()
{
    local shell=$1
    $shell bash_tests/exec14.sub
}

test_quotearray4.sub()
{
    local shell=$1
    $shell bash_tests/quotearray4.sub
}

test_array16.sub()
{
    local shell=$1
    $shell bash_tests/array16.sub
}

test_array17.sub()
{
    local shell=$1
    $shell bash_tests/array17.sub
}

test_quotearray5.sub()
{
    local shell=$1
    $shell bash_tests/quotearray5.sub
}

test_new-exp1.sub()
{
    local shell=$1
    $shell bash_tests/new-exp1.sub
}

test_mapfile2.sub()
{
    local shell=$1
    $shell bash_tests/mapfile2.sub
}

test_dollar-at-star5.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star5.sub
}

test_varenv8.sub()
{
    local shell=$1
    $shell bash_tests/varenv8.sub
}

test_source6.sub()
{
    local shell=$1
    $shell bash_tests/source6.sub
}

test_getopts5.sub()
{
    local shell=$1
    $shell bash_tests/getopts5.sub
}

test_errors7.sub()
{
    local shell=$1
    $shell bash_tests/errors7.sub
}

test_nameref19.sub()
{
    local shell=$1
    $shell bash_tests/nameref19.sub
}

test_printf2.sub()
{
    local shell=$1
    $shell bash_tests/printf2.sub
}

test_nquote2.tests()
{
    local shell=$1
    $shell bash_tests/nquote2.tests
}

test_posixexp1.sub()
{
    local shell=$1
    $shell bash_tests/posixexp1.sub
}

test_exec6.sub()
{
    local shell=$1
    $shell bash_tests/exec6.sub
}

test_exec4.sub()
{
    local shell=$1
    $shell bash_tests/exec4.sub
}

test_posixexp3.sub()
{
    local shell=$1
    $shell bash_tests/posixexp3.sub
}

test_errors5.sub()
{
    local shell=$1
    $shell bash_tests/errors5.sub
}

test_ifs1.sub()
{
    local shell=$1
    $shell bash_tests/ifs1.sub
}

test_getopts7.sub()
{
    local shell=$1
    $shell bash_tests/getopts7.sub
}

test_coproc.tests()
{
    local shell=$1
    $shell bash_tests/coproc.tests
}

test_source4.sub()
{
    local shell=$1
    $shell bash_tests/source4.sub
}

test_new-exp3.sub()
{
    local shell=$1
    $shell bash_tests/new-exp3.sub
}

test_dollar-at-star7.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star7.sub
}

test_array29.sub()
{
    local shell=$1
    $shell bash_tests/array29.sub
}

test_array15.sub()
{
    local shell=$1
    $shell bash_tests/array15.sub
}

test_exportfunc.tests()
{
    local shell=$1
    $shell bash_tests/exportfunc.tests
}

test_array14.sub()
{
    local shell=$1
    $shell bash_tests/array14.sub
}

test_array28.sub()
{
    local shell=$1
    $shell bash_tests/array28.sub
}

test_dollar-at-star6.sub()
{
    local shell=$1
    $shell bash_tests/dollar-at-star6.sub
}

test_new-exp2.sub()
{
    local shell=$1
    $shell bash_tests/new-exp2.sub
}

test_mapfile1.sub()
{
    local shell=$1
    $shell bash_tests/mapfile1.sub
}

test_source5.sub()
{
    local shell=$1
    $shell bash_tests/source5.sub
}

test_getopts6.sub()
{
    local shell=$1
    $shell bash_tests/getopts6.sub
}

test_errors4.sub()
{
    local shell=$1
    $shell bash_tests/errors4.sub
}

test_printf1.sub()
{
    local shell=$1
    $shell bash_tests/printf1.sub
}

test_posixexp2.sub()
{
    local shell=$1
    $shell bash_tests/posixexp2.sub
}

test_exec5.sub()
{
    local shell=$1
    $shell bash_tests/exec5.sub
}

test_heredoc.tests()
{
    local shell=$1
    $shell bash_tests/heredoc.tests
}

test_select.sh()
{
    local shell=$1
    echo 1 | $shell bash_tests/select.sh
}

## We run all tests composed with && to exit on the first that fails
# commented out IFS tests should work but this is an issue with PaSh, bug report made
if [ "$#" -eq 0 ] || [ "$test_mode" = "bash" ]; then
    if [ "$#" -eq 0 ]; then 
        echo "Warning: these tests should be run with bash as the first argument to test the bash mode."
    fi
    # run_test test_exec8.sub - any uncommented script here and forward has an alias (or declare -A) or unset in it
    # run_test test_arith-for.tests - error script, looks good besides an extra error print, this is an important test
    # run_test test_history2.sub - set -o history won't work here bc the history is gonna look different
    # run_test test_posixexp2.tests - setting posix mode won't work this affects parsing
    run_test test_new-exp12.sub 
    run_test test_array4.sub
    run_test test_varenv18.sub
    run_test test_assoc14.sub
    run_test test_attr.tests 
    run_test test_array28.sub
    run_test test_heredoc.tests
    run_test test_jobs1.sub
    run_test test_case1.sub
    # run_test test_nameref17.sub
    # run_test test_set-e3a.sub - echo $- shouldn't work because it depends on variable state
    run_test test_lastpipe1.sub
    run_test test_dollar-at7.sub
    run_test test_trap4.sub
    # run_test test_errors9.sub - this is fine because we're just printing error statements?, dash errors here
    # run_test test_dollar-star3.sub - setting IFS, also doesn't work in dash
    run_test test_assoc1.sub
    # run_test test_alias4.sub
    # run_test test_varenv6.sub
    run_test test_builtins2.sub
    # run_test test_array25.sub
    run_test test_nquote1.sub
    # run_test test_array19.sub
    # run_test test_array18.sub
    run_test test_redir1.sub
    # run_test test_array24.sub
    # run_test test_array30.sub
    run_test test_builtins3.sub
    # run_test test_varenv7.sub
    # run_test test_alias5.sub
    # run_test test_dollar-star2.sub - setting IFS, also doesn't work in dash
    # run_test test_nameref1.sub
    run_test test_errors8.sub
    run_test test_trap5.sub
    run_test test_dollar-at6.sub
    # run_test test_nameref16.sub
    run_test test_exec9.sub
    run_test test_heredoc1.sub
    # run_test test_read1.sub
    run_test test_comsub-eof.tests
    run_test test_jobs2.sub
    run_test test_read3.sub
    # run_test test_heredoc3.sub - error script, pash fails on parsing while bash fails live as expected
    run_test test_trap2a.sub
    run_test test_case2.sub
    # run_test test_nameref14.sub
    # run_test test_lastpipe2.sub
    run_test test_dollar-at4.sub
    # run_test test_nameref3.sub
    # run_test test_quotearray.tests
    run_test test_getopts8.sub
    # run_test test_assoc2.sub
    # run_test test_builtins1.sub
    # run_test test_dollar-at-star8.sub
    # run_test test_varenv5.sub
    # run_test test_type4.sub
    # run_test test_array26.sub
    run_test test_nquote2.sub
    run_test test_redir3.sub
    run_test test_redir2.sub
    run_test test_nquote3.sub
    # run_test test_array27.sub
    # run_test test_varenv4.sub
    # run_test test_dollar-at-star9.sub
    # run_test test_assoc3.sub
    # run_test test_alias6.sub
    # run_test test_posix2syntax.sub
    # run_test test_dollar-star1.sub - setting IFS, also doesn't work in dash
    run_test test_getopts9.sub
    run_test test_dstack2.tests
    # run_test test_trap6.sub
    run_test test_nameref2.sub
    # run_test test_dollar-at5.sub
    # run_test test_arith.tests
    # run_test test_lastpipe3.sub
    # run_test test_nameref15.sub
    run_test test_case3.sub
    # run_test test_procsub.tests
    run_test test_heredoc2.sub
    # run_test test_read2.sub
    run_test test_jobs3.sub
    run_test test_heredoc6.sub
    run_test test_jobs7.sub
    run_test test_read6.sub
    # run_test test_shopt.tests
    # run_test test_comsub.tests
    # run_test test_dollar-star10.sub
    run_test test_dollar-at1.sub
    # run_test test_appendop1.sub - uses a typeset -A
    # run_test test_nameref11.sub
    run_test test_strip.tests
    # run_test test_nameref6.sub
    # run_test test_trap2.sub
    # run_test test_alias2.sub
    run_test test_assoc7.sub
    # run_test test_dollar-star5.sub - setting IFS, also doesn't work in dash
    # run_test test_builtins4.sub
    # run_test test_new-exp9.sub - local -A
    run_test test_redir6.sub
    # run_test test_array23.sub - error script
    run_test test_type1.sub
    # run_test test_array22.sub
    # run_test test_redir7.sub
    run_test test_intl1.sub
    run_test test_varenv1.sub
    # run_test test_builtins5.sub
    # run_test test_new-exp8.sub
    run_test test_dollar-star4.sub
    # run_test test_assoc6.sub
    # run_test test_alias3.sub
    # run_test test_set-x1.sub
    # run_test test_trap3.sub - weird stuff with PS4, also fails without bash mode
    # run_test test_nameref7.sub
    # run_test test_history.tests
    # run_test test_nameref10.sub
    run_test test_parser1.sub
    # run_test test_jobs.tests - error script
    run_test test_posixexp8.sub
    # run_test test_read7.sub
    run_test test_jobs6.sub
    run_test test_heredoc7.sub
    run_test test_heredoc5.sub
    run_test test_read5.sub
    # run_test test_jobs4.sub - also fails with dash pash, jobs get tricky
    # run_test test_case4.sub
    run_test test_dollar-at2.sub
    # run_test test_nameref.tests
    # run_test test_nameref12.sub
    # run_test test_appendop2.sub
    run_test test_trap1.sub
    # run_test test_nameref5.sub
    # run_test test_assoc4.sub - error script
    # run_test test_tilde2.tests
    run_test test_vredir.tests
    # run_test test_alias1.sub
    # run_test test_mapfile.tests
    run_test test_dollar-star6.sub
    run_test test_parser.tests
    # run_test test_varenv3.sub
    run_test test_builtins7.sub
    # run_test test_intl3.sub - error script
    run_test test_redir5.sub
    # run_test test_errors.tests
    run_test test_type2.sub
    # run_test test_array20.sub
    # run_test test_nquote4.sub
    run_test test_type3.sub
    run_test test_nquote5.sub
    # run_test test_array21.sub
    # run_test test_redir4.sub
    # run_test test_intl2.sub
    # run_test test_herestr1.sub - messing with IFS, never good
    # run_test test_builtins6.sub
    run_test test_varenv2.sub
    # run_test test_dollar-star7.sub
    run_test test_assoc5.sub
    # run_test test_nameref4.sub
    # run_test test_nameref13.sub
    run_test test_dollar-at3.sub
    run_test test_extglob3.tests
    # run_test test_jobs5.sub
    run_test test_read4.sub
    run_test test_heredoc4.sub
    # run_test test_dbg-support3.sub - using BASH_ARGV won't work bc call stacks will be different
    run_test test_globstar2.sub
    run_test test_arith4.sub
    run_test test_vredir4.sub
    run_test test_exp6.sub
    run_test test_comsub-eof3.sub
    # run_test test_varenv.tests
    run_test test_rhs-exp.tests
    # run_test test_varenv12.sub
    # run_test test_attr1.sub - alias stuff
    run_test test_glob1.sub
    # run_test test_glob.tests
    # run_test test_varenv13.sub
    run_test test_comsub-eof2.sub
    # run_test test_exp7.sub
    # run_test test_braces.tests
    run_test test_vredir5.sub
    run_test test_globstar3.sub
    run_test test_arith5.sub
    run_test test_vredir7.sub
    run_test test_arith7.sub
    run_test test_globstar1.sub
    run_test test_comsub-eof0.sub
    run_test test_exp5.sub
    # run_test test_varenv11.sub
    # run_test test_glob3.sub
    run_test test_glob10.sub
    # run_test test_dbg-support2.tests - again stuff gets weird when using these BASH variables
    # run_test test_quote4.sub
    # run_test test_dynvar.tests
    # run_test test_case.tests
    run_test test_attr2.sub
    # run_test test_extglob1a.sub - calls shopt -s extglob which changes the behavior of the shell
    # run_test test_glob2.sub
    # run_test test_varenv10.sub
    run_test test_exp4.sub
    run_test test_comsub-eof1.sub
    # run_test test_arith6.sub
    # run_test test_vredir6.sub
    # run_test test_assoc18.sub
    # run_test test_shopt1.sub - again weird shopt stuff here
    # run_test test_array8.sub - error script
    run_test test_vredir2.sub
    run_test test_arith2.sub
    # run_test test_varenv14.sub
    run_test test_comsub-eof5.sub
    # run_test test_rsh1.sub - error script
    run_test test_glob6.sub
    # run_test test_quote1.sub
    # run_test test_read.tests
    # run_test test_cond-regexp1.sub
    run_test test_extglob2.tests
    # run_test test_appendop.tests
    run_test test_iquote1.sub
    run_test test_exportfunc1.sub
    # run_test test_extglob.tests - has shopt -s extglob
    run_test test_glob7.sub
    # run_test test_dbg-support.tests - using BASH variables related to debugging
    run_test test_exp1.sub
    run_test test_comsub-eof4.sub
    # run_test test_rsh.tests - error script
    run_test test_varenv15.sub
    # run_test test_func.tests
    # run_test test_arith3.sub
    # run_test test_vredir3.sub
    # run_test test_rhs-exp1.sub
    # run_test test_array9.sub
    # run_test test_posixexp.tests
    run_test test_arith1.sub
    run_test test_vredir1.sub
    # run_test test_cond.tests
    run_test test_ifs.tests
    # run_test test_varenv17.sub
    # run_test test_exp3.sub - messing with IFS
    run_test test_comsub-eof6.sub
    # run_test test_rsh2.sub - error script
    # run_test test_glob5.sub
    run_test test_quote2.sub
    # run_test test_lastpipe.tests
    # run_test test_exportfunc3.sub
    run_test test_cond-regexp2.sub
    # run_test test_more-exp.tests
    run_test test_cond-regexp3.sub
    run_test test_exportfunc2.sub
    run_test test_quote3.sub
    # run_test test_glob4.sub
    # run_test test_exp2.sub
    # run_test test_varenv16.sub
    # run_test test_assoc13.sub
    run_test test_array7.sub
    # run_test test_exp.tests
    # run_test test_new-exp.tests
    run_test test_exp12.sub
    # run_test test_comsub6.sub
    # run_test test_histexp6.sub
    run_test test_getopts.tests
    run_test test_glob9.sub
    run_test test_input-line.sub
    # run_test test_history3.sub
    run_test test_func1.sub
    # run_test test_extglob4.sub - extend globbing
    # run_test test_dstack.tests
    run_test test_extglob5.sub
    # run_test test_unicode1.sub
    run_test test_glob8.sub
    run_test test_histexp7.sub
    # run_test test_exp13.sub
    run_test test_new-exp13.sub
    # run_test test_assoc.tests
    # run_test test_dbg-support.sub - uses BASH variables doesn't work
    # run_test test_array6.sub
    # run_test test_assoc12.sub
    run_test test_assoc10.sub
    run_test test_posixpipe.tests
    run_test test_new-exp11.sub
    # run_test test_exp11.sub
    # run_test test_comsub5.sub
    # run_test test_posix2.tests
    # run_test test_histexp5.sub - history expansion issue
    # run_test test_printf.tests
    run_test test_unicode3.sub
    run_test test_set-x.tests
    # run_test test_extglob7.sub - extended globbing again
    run_test test_func2.sub
    # run_test test_comsub-posix.tests
    run_test test_quote.tests
    # run_test test_extglob6.sub - extended globbing again
    run_test test_func3.sub
    # run_test test_history1.sub
    run_test test_unicode2.sub
    run_test test_tilde3.sub
    run_test test_casemod.tests
    # run_test test_histexp4.sub - history expansion
    # run_test test_comsub4.sub
    # run_test test_getopts10.sub
    # run_test test_varenv19.sub - same local function issue
    # run_test test_exp10.sub - setting IFS doesn't work
    # run_test test_nquote5.tests
    # run_test test_new-exp10.sub
    # run_test test_array5.sub
    run_test test_test1.sub
    # run_test test_assoc11.sub
    run_test test_array1.sub
    # run_test test_assoc15.sub
    run_test test_new-exp14.sub
    # run_test test_varenv21.sub
    # run_test test_exp9.sub
    run_test test_procsub2.sub
    # run_test test_histexp.tests
    # run_test test_history5.sub
    # run_test test_nquote.tests
    run_test test_extglob2.sub
    run_test test_nquote1.tests
    run_test test_set-e1.sub
    # run_test test_extglob3.sub - extended globbing
    # run_test test_globstar.tests
    run_test test_history4.sub
    # run_test test_intl.tests
    run_test test_comsub1.sub
    run_test test_histexp1.sub
    # run_test test_varenv20.sub
    # run_test test_exp8.sub
    run_test test_new-exp15.sub
    # run_test test_complete.tests
    run_test test_array2.sub
    run_test test_assoc16.sub
    run_test test_vredir8.sub
    run_test test_cond-regexp.sub
    # run_test test_arith8.sub
    run_test test_set-e.tests
    # run_test test_histexp3.sub - history expansion
    # run_test test_varenv22.sub - error script
    run_test test_procsub1.sub
    run_test test_comsub3.sub
    # run_test test_history6.sub
    # run_test test_builtins.tests
    # run_test test_extglob1.sub - extended globbing shouldn't work
    # run_test test_func4.sub
    # run_test test_alias.tests
    run_test test_set-e3.sub
    run_test test_set-e2.sub
    run_test test_nquote3.tests
    run_test test_comsub2.sub
    # run_test test_histexp2.sub - history expansion shouldn't work
    run_test test_new-exp16.sub
    # run_test test_assoc17.sub
    run_test test_array3.sub
    # run_test test_type.tests
    run_test test_exec1.sub
    # run_test test_posixexp6.sub
    # run_test test_redir11.sub
    # run_test test_nameref22.sub
    run_test test_getopts2.sub
    run_test test_nameref9.sub
    # run_test test_assoc8.sub - declare -A here
    run_test test_new-exp6.sub
    # run_test test_dollar-at-star2.sub - IFS stuff again
    run_test test_exec12.sub
    run_test test_source1.sub
    run_test test_array10.sub
    run_test test_redir9.sub
    # run_test test_quotearray2.sub
    run_test test_redir8.sub
    # run_test test_quotearray3.sub
    # run_test test_array11.sub
    # run_test test_dollar-at-star3.sub
    run_test test_exec13.sub
    run_test test_new-exp7.sub
    # run_test test_tilde.tests
    # run_test test_assoc9.sub
    run_test test_nquote4.tests
    run_test test_errors1.sub
    # run_test test_nameref8.sub
    run_test test_getopts3.sub
    run_test test_ifs-posix.tests
    # run_test test_comsub-posix1.sub - error script
    # run_test test_nameref23.sub
    run_test test_redir10.sub
    run_test test_printf4.sub
    run_test test_posixexp7.sub
    run_test test_read8.sub
    run_test test_exec2.sub
    # run_test test_posixexp5.sub
    # run_test test_trap.tests
    # run_test test_redir.tests - error script
    run_test test_redir12.sub
    # run_test test_nameref21.sub
    run_test test_comsub-posix3.sub
    run_test test_getopts1.sub
    # run_test test_dollar-at-star10.sub
    run_test test_errors3.sub
    # run_test test_dollar-star9.sub
    run_test test_exec11.sub
    run_test test_dollar-at-star1.sub
    # run_test test_new-exp5.sub - messes with IFS
    run_test test_source2.sub
    # run_test test_quotearray1.sub
    # run_test test_array13.sub - declare -a
    run_test test_posixpat.tests
    # run_test test_array.tests
    # run_test test_array12.sub
    # run_test test_iquote.tests
    run_test test_source3.sub
    run_test test_new-exp4.sub
    run_test test_cprint.tests
    # run_test test_exec10.sub - error script
    # run_test test_test.tests
    # run_test test_dollar-star8.sub
    run_test test_errors2.sub
    # run_test test_dollar-at-star11.sub
    run_test test_comsub-posix2.sub
    # run_test test_nameref20.sub
    run_test test_invert.tests
    # run_test test_posixexp4.sub
    run_test test_exec3.sub
    run_test test_exec7.sub
    run_test test_precedence.tests
    # run_test test_printf3.sub
    # run_test test_nameref18.sub
    run_test test_comsub-posix6.sub
    # run_test test_errors6.sub
    run_test test_getopts4.sub
    # run_test test_herestr.tests
    # run_test test_source7.sub
    # run_test test_varenv9.sub
    # run_test test_dollar-at-star4.sub
    run_test test_exec14.sub
    # run_test test_quotearray4.sub
    # run_test test_array16.sub
    # run_test test_array17.sub
    # run_test test_quotearray5.sub
    run_test test_new-exp1.sub
    run_test test_mapfile2.sub
    run_test test_dollar-at-star5.sub
    # run_test test_varenv8.sub
    run_test test_source6.sub
    run_test test_getopts5.sub
    run_test test_errors7.sub
    # run_test test_nameref19.sub
    # run_test test_printf2.sub
    run_test test_nquote2.tests
    # run_test test_posixexp1.sub
    run_test test_exec6.sub
    run_test test_exec4.sub
    # run_test test_posixexp3.sub
    run_test test_errors5.sub
    # run_test test_ifs1.sub - messes with IFS
    run_test test_getopts7.sub
    # run_test test_coproc.tests - error script
    run_test test_source4.sub
    run_test test_new-exp3.sub
    # run_test test_dollar-at-star7.sub
    # run_test test_array29.sub
    # run_test test_array15.sub
    # run_test test_exportfunc.tests
    # run_test test_array14.sub
    # run_test test_dollar-at-star6.sub - setting IFS doesn't work
    run_test test_new-exp2.sub
    run_test test_mapfile1.sub
    # run_test test_source5.sub
    run_test test_getopts6.sub
    run_test test_errors4.sub
    # run_test test_printf1.sub
    # run_test test_posixexp2.sub
    run_test test_exec5.sub
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
