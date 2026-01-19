#!/bin/bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
# time: print real in seconds, to simplify parsing

bash="bash"
pash="$PASH_TOP/pa.sh --assert_compiler_success --parallel_pipelines --profile_driven --bash"

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
    $shell jobs1.sub
}

test_exec8.sub()
{
    local shell=$1
    $shell exec8.sub
}

test_case1.sub()
{
    local shell=$1
    $shell case1.sub
}

test_nameref17.sub()
{
    local shell=$1
    $shell nameref17.sub
}

test_lastpipe1.sub()
{
    local shell=$1
    $shell lastpipe1.sub
}

test_dollar-at7.sub()
{
    local shell=$1
    $shell dollar-at7.sub
}

test_trap4.sub()
{
    local shell=$1
    $shell trap4.sub
}

test_errors9.sub()
{
    local shell=$1
    $shell errors9.sub
}

test_dollar-star3.sub()
{
    local shell=$1
    $shell dollar-star3.sub
}

test_assoc1.sub()
{
    local shell=$1
    $shell assoc1.sub
}

test_alias4.sub()
{
    local shell=$1
    $shell alias4.sub
}

test_varenv6.sub()
{
    local shell=$1
    $shell varenv6.sub
}

test_builtins2.sub()
{
    local shell=$1
    $shell builtins2.sub
}

test_array25.sub()
{
    local shell=$1
    $shell array25.sub
}

test_nquote1.sub()
{
    local shell=$1
    $shell nquote1.sub
}

test_array19.sub()
{
    local shell=$1
    $shell array19.sub
}

test_array18.sub()
{
    local shell=$1
    $shell array18.sub
}

test_redir1.sub()
{
    local shell=$1
    $shell redir1.sub
}

test_array24.sub()
{
    local shell=$1
    $shell array24.sub
}

test_array30.sub()
{
    local shell=$1
    $shell array30.sub
}

test_builtins3.sub()
{
    local shell=$1
    $shell builtins3.sub
}

test_varenv7.sub()
{
    local shell=$1
    $shell varenv7.sub
}

test_alias5.sub()
{
    local shell=$1
    $shell alias5.sub
}

test_dollar-star2.sub()
{
    local shell=$1
    $shell dollar-star2.sub
}

test_nameref1.sub()
{
    local shell=$1
    $shell nameref1.sub
}

test_errors8.sub()
{
    local shell=$1
    $shell errors8.sub
}

test_trap5.sub()
{
    local shell=$1
    $shell trap5.sub
}

test_dollar-at6.sub()
{
    local shell=$1
    $shell dollar-at6.sub
}

test_nameref16.sub()
{
    local shell=$1
    $shell nameref16.sub
}

test_exec9.sub()
{
    local shell=$1
    $shell exec9.sub
}

test_heredoc1.sub()
{
    local shell=$1
    $shell heredoc1.sub
}

test_read1.sub()
{
    local shell=$1
    $shell read1.sub
}

test_comsub-eof.tests()
{
    local shell=$1
    $shell comsub-eof.tests
}

test_jobs2.sub()
{
    local shell=$1
    $shell jobs2.sub
}

test_read3.sub()
{
    local shell=$1
    $shell read3.sub
}

test_heredoc3.sub()
{
    local shell=$1
    $shell heredoc3.sub
}

test_trap2a.sub()
{
    local shell=$1
    $shell trap2a.sub
}

test_case2.sub()
{
    local shell=$1
    $shell case2.sub
}

test_nameref14.sub()
{
    local shell=$1
    $shell nameref14.sub
}

test_lastpipe2.sub()
{
    local shell=$1
    $shell lastpipe2.sub
}

test_dollar-at4.sub()
{
    local shell=$1
    $shell dollar-at4.sub
}

test_nameref3.sub()
{
    local shell=$1
    $shell nameref3.sub
}

test_quotearray.tests()
{
    local shell=$1
    $shell quotearray.tests
}

test_getopts8.sub()
{
    local shell=$1
    $shell getopts8.sub
}

test_assoc2.sub()
{
    local shell=$1
    $shell assoc2.sub
}

test_builtins1.sub()
{
    local shell=$1
    $shell builtins1.sub
}

test_dollar-at-star8.sub()
{
    local shell=$1
    $shell dollar-at-star8.sub
}

test_varenv5.sub()
{
    local shell=$1
    $shell varenv5.sub
}

test_type4.sub()
{
    local shell=$1
    $shell type4.sub
}

test_array26.sub()
{
    local shell=$1
    $shell array26.sub
}

test_nquote2.sub()
{
    local shell=$1
    $shell nquote2.sub
}

test_redir3.sub()
{
    local shell=$1
    echo hi1 hi2 hi3 hi4 | $shell redir3.sub
}

test_redir2.sub()
{
    local shell=$1
    $shell redir2.sub
}

test_nquote3.sub()
{
    local shell=$1
    $shell nquote3.sub
}

test_array27.sub()
{
    local shell=$1
    $shell array27.sub
}

test_varenv4.sub()
{
    local shell=$1
    $shell varenv4.sub
}

test_dollar-at-star9.sub()
{
    local shell=$1
    $shell dollar-at-star9.sub
}

test_assoc3.sub()
{
    local shell=$1
    $shell assoc3.sub
}

test_alias6.sub()
{
    local shell=$1
    $shell alias6.sub
}

test_posix2syntax.sub()
{
    local shell=$1
    $shell posix2syntax.sub
}

test_dollar-star1.sub()
{
    local shell=$1
    $shell dollar-star1.sub
}

test_getopts9.sub()
{
    local shell=$1
    $shell getopts9.sub
}

test_dstack2.tests()
{
    local shell=$1
    $shell dstack2.tests
}

test_trap6.sub()
{
    local shell=$1
    $shell trap6.sub
}

test_nameref2.sub()
{
    local shell=$1
    $shell nameref2.sub
}

test_dollar-at5.sub()
{
    local shell=$1
    $shell dollar-at5.sub
}

test_arith.tests()
{
    local shell=$1
    $shell arith.tests
}

test_lastpipe3.sub()
{
    local shell=$1
    $shell lastpipe3.sub
}

test_nameref15.sub()
{
    local shell=$1
    $shell nameref15.sub
}

test_case3.sub()
{
    local shell=$1
    $shell case3.sub
}

test_procsub.tests()
{
    local shell=$1
    $shell procsub.tests
}

test_heredoc2.sub()
{
    local shell=$1
    $shell heredoc2.sub
}

test_read2.sub()
{
    local shell=$1
    $shell read2.sub
}

test_jobs3.sub()
{
    local shell=$1
    $shell jobs3.sub
}

test_heredoc6.sub()
{
    local shell=$1
    $shell heredoc6.sub
}

test_jobs7.sub()
{
    local shell=$1
    $shell jobs7.sub
}

test_read6.sub()
{
    local shell=$1
    $shell read6.sub
}

test_shopt.tests()
{
    local shell=$1
    $shell shopt.tests
}

test_comsub.tests()
{
    local shell=$1
    $shell comsub.tests
}

test_dollar-star10.sub()
{
    local shell=$1
    $shell dollar-star10.sub
}

test_dollar-at1.sub()
{
    local shell=$1
    $shell dollar-at1.sub
}

test_appendop1.sub()
{
    local shell=$1
    $shell appendop1.sub
}

test_nameref11.sub()
{
    local shell=$1
    $shell nameref11.sub
}

test_strip.tests()
{
    local shell=$1
    $shell strip.tests
}

test_nameref6.sub()
{
    local shell=$1
    $shell nameref6.sub
}

test_trap2.sub()
{
    local shell=$1
    $shell trap2.sub
}

test_alias2.sub()
{
    local shell=$1
    $shell alias2.sub
}

test_assoc7.sub()
{
    local shell=$1
    $shell assoc7.sub
}

test_dollar-star5.sub()
{
    local shell=$1
    $shell dollar-star5.sub
}

test_builtins4.sub()
{
    local shell=$1
    $shell builtins4.sub
}

test_new-exp9.sub()
{
    local shell=$1
    $shell new-exp9.sub
}

test_redir6.sub()
{
    local shell=$1
    $shell redir6.sub
}

test_array23.sub()
{
    local shell=$1
    $shell array23.sub
}

test_type1.sub()
{
    local shell=$1
    $shell type1.sub
}

test_array22.sub()
{
    local shell=$1
    $shell array22.sub
}

test_redir7.sub()
{
    local shell=$1
    $shell redir7.sub
}

test_intl1.sub()
{
    local shell=$1
    $shell intl1.sub
}

test_varenv1.sub()
{
    local shell=$1
    $shell varenv1.sub
}

test_builtins5.sub()
{
    local shell=$1
    $shell builtins5.sub
}

test_new-exp8.sub()
{
    local shell=$1
    $shell new-exp8.sub
}

test_dollar-star4.sub()
{
    local shell=$1
    $shell dollar-star4.sub
}

test_assoc6.sub()
{
    local shell=$1
    $shell assoc6.sub
}

test_alias3.sub()
{
    local shell=$1
    $shell alias3.sub
}

test_set-x1.sub()
{
    local shell=$1
    $shell set-x1.sub
}

test_trap3.sub()
{
    local shell=$1
    $shell trap3.sub
}

test_nameref7.sub()
{
    local shell=$1
    $shell nameref7.sub
}

test_history.tests()
{
    local shell=$1
    $shell history.tests
}

test_nameref10.sub()
{
    local shell=$1
    $shell nameref10.sub
}

test_parser1.sub()
{
    local shell=$1
    $shell parser1.sub
}

test_jobs.tests()
{
    local shell=$1
    $shell jobs.tests
}

test_posixexp8.sub()
{
    local shell=$1
    $shell posixexp8.sub
}

test_read7.sub()
{
    local shell=$1
    $shell read7.sub
}

test_jobs6.sub()
{
    local shell=$1
    $shell jobs6.sub
}

test_heredoc7.sub()
{
    local shell=$1
    $shell heredoc7.sub
}

test_heredoc5.sub()
{
    local shell=$1
    $shell heredoc5.sub
}

test_read5.sub()
{
    local shell=$1
    $shell read5.sub
}

test_jobs4.sub()
{
    local shell=$1
    $shell jobs4.sub
}

test_case4.sub()
{
    local shell=$1
    $shell case4.sub
}

test_dollar-at2.sub()
{
    local shell=$1
    $shell dollar-at2.sub
}

test_nameref.tests()
{
    local shell=$1
    $shell nameref.tests
}

test_nameref12.sub()
{
    local shell=$1
    $shell nameref12.sub
}

test_appendop2.sub()
{
    local shell=$1
    $shell appendop2.sub
}

test_trap1.sub()
{
    local shell=$1
    $shell trap1.sub
}

test_nameref5.sub()
{
    local shell=$1
    $shell nameref5.sub
}

test_assoc4.sub()
{
    local shell=$1
    $shell assoc4.sub
}

test_tilde2.tests()
{
    local shell=$1
    $shell tilde2.tests
}

test_vredir.tests()
{
    local shell=$1
    $shell vredir.tests
}

test_alias1.sub()
{
    local shell=$1
    $shell alias1.sub
}

test_mapfile.tests()
{
    local shell=$1
    $shell mapfile.tests
}

test_dollar-star6.sub()
{
    local shell=$1
    $shell dollar-star6.sub
}

test_parser.tests()
{
    local shell=$1
    $shell parser.tests
}

test_varenv3.sub()
{
    local shell=$1
    $shell varenv3.sub
}

test_builtins7.sub()
{
    local shell=$1
    $shell builtins7.sub
}

test_intl3.sub()
{
    local shell=$1
    $shell intl3.sub
}

test_redir5.sub()
{
    local shell=$1
    $shell redir5.sub
}

test_errors.tests()
{
    local shell=$1
    $shell errors.tests
}

test_type2.sub()
{
    local shell=$1
    $shell type2.sub
}

test_array20.sub()
{
    local shell=$1
    $shell array20.sub
}

test_nquote4.sub()
{
    local shell=$1
    $shell nquote4.sub
}

test_type3.sub()
{
    local shell=$1
    $shell type3.sub
}

test_nquote5.sub()
{
    local shell=$1
    $shell nquote5.sub
}

test_array21.sub()
{
    local shell=$1
    $shell array21.sub
}

test_redir4.sub()
{
    local shell=$1
    $shell redir4.sub
}

test_intl2.sub()
{
    local shell=$1
    $shell intl2.sub
}

test_herestr1.sub()
{
    local shell=$1
    $shell herestr1.sub
}

test_builtins6.sub()
{
    local shell=$1
    $shell builtins6.sub
}

test_varenv2.sub()
{
    local shell=$1
    $shell varenv2.sub
}

test_dollar-star7.sub()
{
    local shell=$1
    $shell dollar-star7.sub
}

test_assoc5.sub()
{
    local shell=$1
    $shell assoc5.sub
}

test_nameref4.sub()
{
    local shell=$1
    $shell nameref4.sub
}

test_nameref13.sub()
{
    local shell=$1
    $shell nameref13.sub
}

test_dollar-at3.sub()
{
    local shell=$1
    $shell dollar-at3.sub
}

test_extglob3.tests()
{
    local shell=$1
    $shell extglob3.tests
}

test_jobs5.sub()
{
    local shell=$1
    $shell jobs5.sub
}

test_read4.sub()
{
    local shell=$1
    $shell read4.sub
}

test_heredoc4.sub()
{
    local shell=$1
    $shell heredoc4.sub
}

test_dbg-support3.sub()
{
    local shell=$1
    $shell dbg-support3.sub
}

test_globstar2.sub()
{
    local shell=$1
    $shell globstar2.sub
}

test_arith4.sub()
{
    local shell=$1
    $shell arith4.sub
}

test_vredir4.sub()
{
    local shell=$1
    $shell vredir4.sub
}

test_exp6.sub()
{
    local shell=$1
    $shell exp6.sub
}

test_comsub-eof3.sub()
{
    local shell=$1
    $shell comsub-eof3.sub
}

test_varenv.tests()
{
    local shell=$1
    $shell varenv.tests
}

test_rhs-exp.tests()
{
    local shell=$1
    $shell rhs-exp.tests
}

test_varenv12.sub()
{
    local shell=$1
    $shell varenv12.sub
}

test_attr1.sub()
{
    local shell=$1
    $shell attr1.sub
}

test_glob1.sub()
{
    local shell=$1
    $shell glob1.sub
}

test_glob.tests()
{
    local shell=$1
    $shell glob.tests
}

test_set-e3a.sub()
{
    local shell=$1
    $shell set-e3a.sub
}

test_varenv13.sub()
{
    local shell=$1
    $shell varenv13.sub
}

test_comsub-eof2.sub()
{
    local shell=$1
    $shell comsub-eof2.sub
}

test_exp7.sub()
{
    local shell=$1
    $shell exp7.sub
}

test_braces.tests()
{
    local shell=$1
    $shell braces.tests
}

test_vredir5.sub()
{
    local shell=$1
    $shell vredir5.sub
}

test_globstar3.sub()
{
    local shell=$1
    $shell globstar3.sub
}

test_arith5.sub()
{
    local shell=$1
    $shell arith5.sub
}

test_vredir7.sub()
{
    local shell=$1
    $shell vredir7.sub
}

test_arith7.sub()
{
    local shell=$1
    $shell arith7.sub
}

test_globstar1.sub()
{
    local shell=$1
    $shell globstar1.sub
}

test_comsub-eof0.sub()
{
    local shell=$1
    $shell comsub-eof0.sub
}

test_exp5.sub()
{
    local shell=$1
    $shell exp5.sub
}

test_varenv11.sub()
{
    local shell=$1
    $shell varenv11.sub
}

test_glob3.sub()
{
    local shell=$1
    $shell glob3.sub
}

test_glob10.sub()
{
    local shell=$1
    $shell glob10.sub
}

test_dbg-support2.tests()
{
    local shell=$1
    $shell dbg-support2.tests
}

test_quote4.sub()
{
    local shell=$1
    $shell quote4.sub
}

test_dynvar.tests()
{
    local shell=$1
    $shell dynvar.tests
}

test_case.tests()
{
    local shell=$1
    $shell case.tests
}

test_attr2.sub()
{
    local shell=$1
    $shell attr2.sub
}

test_extglob1a.sub()
{
    local shell=$1
    $shell extglob1a.sub
}

test_glob2.sub()
{
    local shell=$1
    $shell glob2.sub
}

test_varenv10.sub()
{
    local shell=$1
    $shell varenv10.sub
}

test_exp4.sub()
{
    local shell=$1
    $shell exp4.sub
}

test_comsub-eof1.sub()
{
    local shell=$1
    $shell comsub-eof1.sub
}

test_arith6.sub()
{
    local shell=$1
    $shell arith6.sub
}

test_vredir6.sub()
{
    local shell=$1
    $shell vredir6.sub
}

test_assoc18.sub()
{
    local shell=$1
    $shell assoc18.sub
}

test_shopt1.sub()
{
    local shell=$1
    $shell shopt1.sub
}

test_array8.sub()
{
    local shell=$1
    $shell array8.sub
}

test_vredir2.sub()
{
    local shell=$1
    $shell vredir2.sub
}

test_arith2.sub()
{
    local shell=$1
    $shell arith2.sub
}

test_varenv14.sub()
{
    local shell=$1
    $shell varenv14.sub
}

test_comsub-eof5.sub()
{
    local shell=$1
    $shell comsub-eof5.sub
}

test_rsh1.sub()
{
    local shell=$1
    $shell rsh1.sub
}

test_glob6.sub()
{
    local shell=$1
    $shell glob6.sub
}

test_arith-for.tests()
{
    local shell=$1
    $shell arith-for.tests
}

test_quote1.sub()
{
    local shell=$1
    $shell quote1.sub
}

test_read.tests()
{
    local shell=$1
    $shell read.tests
}

test_cond-regexp1.sub()
{
    local shell=$1
    $shell cond-regexp1.sub
}

test_extglob2.tests()
{
    local shell=$1
    $shell extglob2.tests
}

test_appendop.tests()
{
    local shell=$1
    $shell appendop.tests
}

test_iquote1.sub()
{
    local shell=$1
    $shell iquote1.sub
}

test_exportfunc1.sub()
{
    local shell=$1
    $shell exportfunc1.sub
}

test_extglob.tests()
{
    local shell=$1
    $shell extglob.tests
}

test_glob7.sub()
{
    local shell=$1
    $shell glob7.sub
}

test_dbg-support.tests()
{
    local shell=$1
    $shell dbg-support.tests
}

test_exp1.sub()
{
    local shell=$1
    $shell exp1.sub
}

test_comsub-eof4.sub()
{
    local shell=$1
    $shell comsub-eof4.sub
}

test_rsh.tests()
{
    local shell=$1
    $shell rsh.tests
}

test_varenv15.sub()
{
    local shell=$1
    $shell varenv15.sub
}

test_func.tests()
{
    local shell=$1
    $shell func.tests
}

test_arith3.sub()
{
    local shell=$1
    $shell arith3.sub
}

test_vredir3.sub()
{
    local shell=$1
    $shell vredir3.sub
}

test_rhs-exp1.sub()
{
    local shell=$1
    $shell rhs-exp1.sub
}

test_array9.sub()
{
    local shell=$1
    $shell array9.sub
}

test_posixexp.tests()
{
    local shell=$1
    $shell posixexp.tests
}

test_arith1.sub()
{
    local shell=$1
    $shell arith1.sub
}

test_vredir1.sub()
{
    local shell=$1
    $shell vredir1.sub
}

test_cond.tests()
{
    local shell=$1
    $shell cond.tests
}

test_ifs.tests()
{
    local shell=$1
    $shell ifs.tests
}

test_varenv17.sub()
{
    local shell=$1
    $shell varenv17.sub
}

test_exp3.sub()
{
    local shell=$1
    $shell exp3.sub
}

test_comsub-eof6.sub()
{
    local shell=$1
    $shell comsub-eof6.sub
}

test_rsh2.sub()
{
    local shell=$1
    $shell rsh2.sub
}

test_glob5.sub()
{
    local shell=$1
    $shell glob5.sub
}

test_quote2.sub()
{
    local shell=$1
    $shell quote2.sub
}

test_lastpipe.tests()
{
    local shell=$1
    $shell lastpipe.tests
}

test_exportfunc3.sub()
{
    local shell=$1
    $shell exportfunc3.sub
}

test_cond-regexp2.sub()
{
    local shell=$1
    $shell cond-regexp2.sub
}

test_more-exp.tests()
{
    local shell=$1
    $shell more-exp.tests
}

test_cond-regexp3.sub()
{
    local shell=$1
    $shell cond-regexp3.sub
}

test_exportfunc2.sub()
{
    local shell=$1
    $shell exportfunc2.sub
}

test_quote3.sub()
{
    local shell=$1
    $shell quote3.sub
}

test_glob4.sub()
{
    local shell=$1
    $shell glob4.sub
}

test_exp2.sub()
{
    local shell=$1
    $shell exp2.sub
}

test_varenv16.sub()
{
    local shell=$1
    $shell varenv16.sub
}

test_assoc13.sub()
{
    local shell=$1
    $shell assoc13.sub
}

test_array7.sub()
{
    local shell=$1
    $shell array7.sub
}

test_exp.tests()
{
    local shell=$1
    $shell exp.tests
}

test_new-exp12.sub()
{
    local shell=$1
    $shell new-exp12.sub
}

test_new-exp.tests()
{
    local shell=$1
    $shell new-exp.tests
}

test_exp12.sub()
{
    local shell=$1
    $shell exp12.sub
}

test_comsub6.sub()
{
    local shell=$1
    $shell comsub6.sub
}

test_histexp6.sub()
{
    local shell=$1
    $shell histexp6.sub
}

test_getopts.tests()
{
    local shell=$1
    $shell getopts.tests
}

test_glob9.sub()
{
    local shell=$1
    $shell glob9.sub
}

test_input-line.sub()
{
    local shell=$1
    echo i | $shell input-line.sub
}

test_history3.sub()
{
    local shell=$1
    $shell history3.sub
}

test_func1.sub()
{
    local shell=$1
    $shell func1.sub
}

test_extglob4.sub()
{
    local shell=$1
    $shell extglob4.sub
}

test_dstack.tests()
{
    local shell=$1
    $shell dstack.tests
}

test_extglob5.sub()
{
    local shell=$1
    $shell extglob5.sub
}

test_history2.sub()
{
    local shell=$1
    $shell history2.sub
}

test_unicode1.sub()
{
    local shell=$1
    $shell unicode1.sub
}

test_glob8.sub()
{
    local shell=$1
    $shell glob8.sub
}

test_histexp7.sub()
{
    local shell=$1
    $shell histexp7.sub
}

test_exp13.sub()
{
    local shell=$1
    $shell exp13.sub
}

test_new-exp13.sub()
{
    local shell=$1
    $shell new-exp13.sub
}

test_assoc.tests()
{
    local shell=$1
    $shell assoc.tests
}

test_dbg-support.sub()
{
    local shell=$1
    $shell dbg-support.sub
}

test_array6.sub()
{
    local shell=$1
    $shell array6.sub
}

test_assoc12.sub()
{
    local shell=$1
    $shell assoc12.sub
}

test_assoc10.sub()
{
    local shell=$1
    $shell assoc10.sub
}

test_array4.sub()
{
    local shell=$1
    $shell array4.sub
}

test_posixpipe.tests()
{
    local shell=$1
    $shell posixpipe.tests
}

test_new-exp11.sub()
{
    local shell=$1
    $shell new-exp11.sub
}

test_exp11.sub()
{
    local shell=$1
    $shell exp11.sub
}

test_comsub5.sub()
{
    local shell=$1
    $shell comsub5.sub
}

test_varenv18.sub()
{
    local shell=$1
    $shell varenv18.sub
}

test_posix2.tests()
{
    local shell=$1
    $shell posix2.tests
}

test_histexp5.sub()
{
    local shell=$1
    $shell histexp5.sub
}

test_printf.tests()
{
    local shell=$1
    $shell printf.tests
}

test_unicode3.sub()
{
    local shell=$1
    $shell unicode3.sub
}

test_set-x.tests()
{
    local shell=$1
    $shell set-x.tests
}

test_extglob7.sub()
{
    local shell=$1
    $shell extglob7.sub
}

test_func2.sub()
{
    local shell=$1
    $shell func2.sub
}

test_comsub-posix.tests()
{
    local shell=$1
    $shell comsub-posix.tests
}

test_quote.tests()
{
    local shell=$1
    $shell quote.tests
}

test_extglob6.sub()
{
    local shell=$1
    $shell extglob6.sub
}

test_func3.sub()
{
    local shell=$1
    $shell func3.sub
}

test_history1.sub()
{
    local shell=$1
    $shell history1.sub
}

test_unicode2.sub()
{
    local shell=$1
    $shell unicode2.sub
}

test_tilde3.sub()
{
    local shell=$1
    $shell tilde3.sub
}

test_casemod.tests()
{
    local shell=$1
    $shell casemod.tests
}

test_histexp4.sub()
{
    local shell=$1
    $shell histexp4.sub
}

test_comsub4.sub()
{
    local shell=$1
    $shell comsub4.sub
}

test_getopts10.sub()
{
    local shell=$1
    $shell getopts10.sub
}

test_varenv19.sub()
{
    local shell=$1
    $shell varenv19.sub
}

test_exp10.sub()
{
    local shell=$1
    $shell exp10.sub
}

test_nquote5.tests()
{
    local shell=$1
    $shell nquote5.tests
}

test_new-exp10.sub()
{
    local shell=$1
    $shell new-exp10.sub
}

test_array5.sub()
{
    local shell=$1
    $shell array5.sub
}

test_test1.sub()
{
    local shell=$1
    $shell test1.sub
}

test_assoc11.sub()
{
    local shell=$1
    $shell assoc11.sub
}

test_array1.sub()
{
    local shell=$1
    $shell array1.sub
}

test_assoc15.sub()
{
    local shell=$1
    $shell assoc15.sub
}

test_new-exp14.sub()
{
    local shell=$1
    $shell new-exp14.sub
}

test_varenv21.sub()
{
    local shell=$1
    $shell varenv21.sub
}

test_exp9.sub()
{
    local shell=$1
    $shell exp9.sub
}

test_procsub2.sub()
{
    local shell=$1
    $shell procsub2.sub
}

test_histexp.tests()
{
    local shell=$1
    $shell histexp.tests
}

test_history5.sub()
{
    local shell=$1
    $shell history5.sub
}

test_nquote.tests()
{
    local shell=$1
    $shell nquote.tests
}

test_extglob2.sub()
{
    local shell=$1
    $shell extglob2.sub
}

test_nquote1.tests()
{
    local shell=$1
    $shell nquote1.tests
}

test_set-e1.sub()
{
    local shell=$1
    $shell set-e1.sub
}

test_extglob3.sub()
{
    local shell=$1
    $shell extglob3.sub
}

test_globstar.tests()
{
    local shell=$1
    $shell globstar.tests
}

test_history4.sub()
{
    local shell=$1
    $shell history4.sub
}

test_intl.tests()
{
    local shell=$1
    $shell intl.tests
}

test_comsub1.sub()
{
    local shell=$1
    $shell comsub1.sub
}

test_histexp1.sub()
{
    local shell=$1
    $shell histexp1.sub
}

test_varenv20.sub()
{
    local shell=$1
    $shell varenv20.sub
}

test_exp8.sub()
{
    local shell=$1
    $shell exp8.sub
}

test_new-exp15.sub()
{
    local shell=$1
    $shell new-exp15.sub
}

test_complete.tests()
{
    local shell=$1
    $shell complete.tests
}

test_assoc14.sub()
{
    local shell=$1
    $shell assoc14.sub
}

test_array2.sub()
{
    local shell=$1
    $shell array2.sub
}

test_assoc16.sub()
{
    local shell=$1
    $shell assoc16.sub
}

test_vredir8.sub()
{
    local shell=$1
    $shell vredir8.sub
}

test_cond-regexp.sub()
{
    local shell=$1
    $shell cond-regexp.sub
}

test_arith8.sub()
{
    local shell=$1
    $shell arith8.sub
}

test_set-e.tests()
{
    local shell=$1
    $shell set-e.tests
}

test_histexp3.sub()
{
    local shell=$1
    $shell histexp3.sub
}

test_varenv22.sub()
{
    local shell=$1
    $shell varenv22.sub
}

test_procsub1.sub()
{
    local shell=$1
    $shell procsub1.sub
}

test_posixexp2.tests()
{
    local shell=$1
    $shell posixexp2.tests
}

test_comsub3.sub()
{
    local shell=$1
    $shell comsub3.sub
}

test_history6.sub()
{
    local shell=$1
    $shell history6.sub
}

test_builtins.tests()
{
    local shell=$1
    $shell builtins.tests
}

test_extglob1.sub()
{
    local shell=$1
    $shell extglob1.sub
}

test_func4.sub()
{
    local shell=$1
    $shell func4.sub
}

test_alias.tests()
{
    local shell=$1
    $shell alias.tests
}

test_set-e3.sub()
{
    local shell=$1
    $shell set-e3.sub
}

test_set-e2.sub()
{
    local shell=$1
    $shell set-e2.sub
}

test_nquote3.tests()
{
    local shell=$1
    $shell nquote3.tests
}

test_comsub2.sub()
{
    local shell=$1
    $shell comsub2.sub
}

test_histexp2.sub()
{
    local shell=$1
    $shell histexp2.sub
}

test_new-exp16.sub()
{
    local shell=$1
    $shell new-exp16.sub
}

test_assoc17.sub()
{
    local shell=$1
    $shell assoc17.sub
}

test_array3.sub()
{
    local shell=$1
    $shell array3.sub
}

test_type.tests()
{
    local shell=$1
    $shell type.tests
}

test_exec1.sub()
{
    local shell=$1
    $shell exec1.sub
}

test_posixexp6.sub()
{
    local shell=$1
    $shell posixexp6.sub
}

test_redir11.sub()
{
    local shell=$1
    $shell redir11.sub
}

test_nameref22.sub()
{
    local shell=$1
    $shell nameref22.sub
}

test_getopts2.sub()
{
    local shell=$1
    $shell getopts2.sub
}

test_nameref9.sub()
{
    local shell=$1
    $shell nameref9.sub
}

test_assoc8.sub()
{
    local shell=$1
    $shell assoc8.sub
}

test_new-exp6.sub()
{
    local shell=$1
    $shell new-exp6.sub
}

test_dollar-at-star2.sub()
{
    local shell=$1
    $shell dollar-at-star2.sub
}

test_exec12.sub()
{
    local shell=$1
    $shell exec12.sub
}

test_source1.sub()
{
    local shell=$1
    $shell source1.sub
}

test_array10.sub()
{
    local shell=$1
    $shell array10.sub
}

test_redir9.sub()
{
    local shell=$1
    $shell redir9.sub
}

test_quotearray2.sub()
{
    local shell=$1
    $shell quotearray2.sub
}

test_redir8.sub()
{
    local shell=$1
    $shell redir8.sub
}

test_quotearray3.sub()
{
    local shell=$1
    $shell quotearray3.sub
}

test_array11.sub()
{
    local shell=$1
    $shell array11.sub
}

test_dollar-at-star3.sub()
{
    local shell=$1
    $shell dollar-at-star3.sub
}

test_exec13.sub()
{
    local shell=$1
    $shell exec13.sub
}

test_new-exp7.sub()
{
    local shell=$1
    $shell new-exp7.sub
}

test_tilde.tests()
{
    local shell=$1
    $shell tilde.tests
}

test_assoc9.sub()
{
    local shell=$1
    $shell assoc9.sub
}

test_nquote4.tests()
{
    local shell=$1
    $shell nquote4.tests
}

test_errors1.sub()
{
    local shell=$1
    $shell errors1.sub
}

test_nameref8.sub()
{
    local shell=$1
    $shell nameref8.sub
}

test_getopts3.sub()
{
    local shell=$1
    $shell getopts3.sub
}

test_ifs-posix.tests()
{
    local shell=$1
    $shell ifs-posix.tests
}

test_comsub-posix1.sub()
{
    local shell=$1
    $shell comsub-posix1.sub
}

test_nameref23.sub()
{
    local shell=$1
    $shell nameref23.sub
}

test_redir10.sub()
{
    local shell=$1
    $shell redir10.sub
}

test_printf4.sub()
{
    local shell=$1
    $shell printf4.sub
}

test_posixexp7.sub()
{
    local shell=$1
    $shell posixexp7.sub
}

test_exec2.sub()
{
    local shell=$1
    $shell exec2.sub
}

test_posixexp5.sub()
{
    local shell=$1
    $shell posixexp5.sub
}

test_trap.tests()
{
    local shell=$1
    $shell trap.tests
}

test_redir.tests()
{
    local shell=$1
    $shell redir.tests
}

test_redir12.sub()
{
    local shell=$1
    $shell redir12.sub
}

test_nameref21.sub()
{
    local shell=$1
    $shell nameref21.sub
}

test_comsub-posix3.sub()
{
    local shell=$1
    $shell comsub-posix3.sub
}

test_getopts1.sub()
{
    local shell=$1
    $shell getopts1.sub
}

test_dollar-at-star10.sub()
{
    local shell=$1
    $shell dollar-at-star10.sub
}

test_errors3.sub()
{
    local shell=$1
    $shell errors3.sub
}

test_dollar-star9.sub()
{
    local shell=$1
    $shell dollar-star9.sub
}

test_exec11.sub()
{
    local shell=$1
    $shell exec11.sub
}

test_dollar-at-star1.sub()
{
    local shell=$1
    $shell dollar-at-star1.sub
}

test_new-exp5.sub()
{
    local shell=$1
    $shell new-exp5.sub
}

test_source2.sub()
{
    local shell=$1
    $shell source2.sub
}

test_quotearray1.sub()
{
    local shell=$1
    $shell quotearray1.sub
}

test_array13.sub()
{
    local shell=$1
    $shell array13.sub
}

test_posixpat.tests()
{
    local shell=$1
    $shell posixpat.tests
}

test_array.tests()
{
    local shell=$1
    $shell array.tests
}

test_array12.sub()
{
    local shell=$1
    $shell array12.sub
}

test_iquote.tests()
{
    local shell=$1
    $shell iquote.tests
}

test_source3.sub()
{
    local shell=$1
    $shell source3.sub
}

test_new-exp4.sub()
{
    local shell=$1
    $shell new-exp4.sub
}

test_cprint.tests()
{
    local shell=$1
    $shell cprint.tests
}

test_exec10.sub()
{
    local shell=$1
    $shell exec10.sub
}

test_test.tests()
{
    local shell=$1
    $shell test.tests
}

test_dollar-star8.sub()
{
    local shell=$1
    $shell dollar-star8.sub
}

test_errors2.sub()
{
    local shell=$1
    $shell errors2.sub
}

test_dollar-at-star11.sub()
{
    local shell=$1
    $shell dollar-at-star11.sub
}

test_comsub-posix2.sub()
{
    local shell=$1
    $shell comsub-posix2.sub
}

test_attr.tests()
{
    local shell=$1
    $shell attr.tests
}

test_nameref20.sub()
{
    local shell=$1
    $shell nameref20.sub
}

test_invert.tests()
{
    local shell=$1
    $shell invert.tests
}

test_posixexp4.sub()
{
    local shell=$1
    $shell posixexp4.sub
}

test_exec3.sub()
{
    local shell=$1
    $shell exec3.sub
}

test_exec7.sub()
{
    local shell=$1
    $shell exec7.sub
}

test_precedence.tests()
{
    local shell=$1
    $shell precedence.tests
}

test_printf3.sub()
{
    local shell=$1
    $shell printf3.sub
}

test_nameref18.sub()
{
    local shell=$1
    $shell nameref18.sub
}

test_comsub-posix6.sub()
{
    local shell=$1
    $shell comsub-posix6.sub
}

test_errors6.sub()
{
    local shell=$1
    $shell errors6.sub
}

test_getopts4.sub()
{
    local shell=$1
    $shell getopts4.sub
}

test_herestr.tests()
{
    local shell=$1
    $shell herestr.tests
}

test_source7.sub()
{
    local shell=$1
    $shell source7.sub
}

test_varenv9.sub()
{
    local shell=$1
    $shell varenv9.sub
}

test_dollar-at-star4.sub()
{
    local shell=$1
    $shell dollar-at-star4.sub
}

test_exec14.sub()
{
    local shell=$1
    $shell exec14.sub
}

test_quotearray4.sub()
{
    local shell=$1
    $shell quotearray4.sub
}

test_array16.sub()
{
    local shell=$1
    $shell array16.sub
}

test_array17.sub()
{
    local shell=$1
    $shell array17.sub
}

test_quotearray5.sub()
{
    local shell=$1
    $shell quotearray5.sub
}

test_new-exp1.sub()
{
    local shell=$1
    $shell new-exp1.sub
}

test_mapfile2.sub()
{
    local shell=$1
    $shell mapfile2.sub
}

test_dollar-at-star5.sub()
{
    local shell=$1
    $shell dollar-at-star5.sub
}

test_varenv8.sub()
{
    local shell=$1
    $shell varenv8.sub
}

test_source6.sub()
{
    local shell=$1
    $shell source6.sub
}

test_getopts5.sub()
{
    local shell=$1
    $shell getopts5.sub
}

test_errors7.sub()
{
    local shell=$1
    $shell errors7.sub
}

test_nameref19.sub()
{
    local shell=$1
    $shell nameref19.sub
}

test_printf2.sub()
{
    local shell=$1
    $shell printf2.sub
}

test_nquote2.tests()
{
    local shell=$1
    $shell nquote2.tests
}

test_posixexp1.sub()
{
    local shell=$1
    $shell posixexp1.sub
}

test_exec6.sub()
{
    local shell=$1
    $shell exec6.sub
}

test_exec4.sub()
{
    local shell=$1
    $shell exec4.sub
}

test_posixexp3.sub()
{
    local shell=$1
    $shell posixexp3.sub
}

test_errors5.sub()
{
    local shell=$1
    $shell errors5.sub
}

test_ifs1.sub()
{
    local shell=$1
    $shell ifs1.sub
}

test_getopts7.sub()
{
    local shell=$1
    $shell getopts7.sub
}

test_coproc.tests()
{
    local shell=$1
    $shell coproc.tests
}

test_source4.sub()
{
    local shell=$1
    $shell source4.sub
}

test_new-exp3.sub()
{
    local shell=$1
    $shell new-exp3.sub
}

test_dollar-at-star7.sub()
{
    local shell=$1
    $shell dollar-at-star7.sub
}

test_array29.sub()
{
    local shell=$1
    $shell array29.sub
}

test_array15.sub()
{
    local shell=$1
    $shell array15.sub
}

test_exportfunc.tests()
{
    local shell=$1
    $shell exportfunc.tests
}

test_array14.sub()
{
    local shell=$1
    $shell array14.sub
}

test_array28.sub()
{
    local shell=$1
    $shell array28.sub
}

test_dollar-at-star6.sub()
{
    local shell=$1
    $shell dollar-at-star6.sub
}

test_new-exp2.sub()
{
    local shell=$1
    $shell new-exp2.sub
}

test_mapfile1.sub()
{
    local shell=$1
    $shell mapfile1.sub
}

test_source5.sub()
{
    local shell=$1
    $shell source5.sub
}

test_getopts6.sub()
{
    local shell=$1
    $shell getopts6.sub
}

test_errors4.sub()
{
    local shell=$1
    $shell errors4.sub
}

test_printf1.sub()
{
    local shell=$1
    $shell printf1.sub
}

test_posixexp2.sub()
{
    local shell=$1
    $shell posixexp2.sub
}

test_exec5.sub()
{
    local shell=$1
    $shell exec5.sub
}

test_heredoc.tests()
{
    local shell=$1
    $shell heredoc.tests
}

test_select.sh()
{
    local shell=$1
    echo 1 | $shell select.sh
}

## We run all tests composed with && to exit on the first that fails
if [ "$#" -eq 0 ] || [ "$test_mode" = "bash" ]; then
    cd "$PASH_TOP/evaluation/tests/interface_tests/bash_tests"
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
    # run_test test_heredoc.tests
    run_test test_jobs1.sub
    run_test test_case1.sub
    # run_test test_nameref17.sub
    # run_test test_set-e3a.sub - echo $- shouldn't work because it depends on variable state
    run_test test_lastpipe1.sub
    run_test test_dollar-at7.sub
    run_test test_trap4.sub
    # run_test test_errors9.sub - this is fine because we're just printing error statements?, dash errors here
    run_test test_dollar-star3.sub
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
    run_test test_dollar-star2.sub
    # run_test test_nameref1.sub
    run_test test_errors8.sub
    # run_test test_trap5.sub - I don't think this should neccesarly work since we could be catching other PaSh signals
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
    run_test test_dollar-star1.sub
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
    # run_test test_heredoc6.sub
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
    run_test test_dollar-star5.sub
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
    # run_test test_posixexp8.sub
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
    # run_test test_vredir.tests - this is passing, however this script assigns file descriptors which will be different
    # depending on what else is open, so output will look different with and without pash
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
    # run_test test_herestr1.sub
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
    # run_test test_vredir4.sub - see test_vredir.tests
    run_test test_exp6.sub
    run_test test_comsub-eof3.sub
    # run_test test_varenv.tests
    run_test test_rhs-exp.tests
    # run_test test_varenv12.sub
    # run_test test_attr1.sub - alias stuff
    run_test test_glob1.sub
    # run_test test_glob.tests
    # run_test test_varenv13.sub
    # run_test test_comsub-eof2.sub
    # run_test test_exp7.sub
    # run_test test_braces.tests
    # run_test test_vredir5.sub - see test_vredir.tests
    run_test test_globstar3.sub
    run_test test_arith5.sub
    # run_test test_vredir7.sub - see test_vredir.tests
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
    # run_test test_vredir2.sub - test_vredir2.sub
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
    # run_test test_exp1.sub -\x7f character causes the expansion server to fail
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
    # run_test test_vredir1.sub - see test_vredir.tests
    # run_test test_cond.tests
    run_test test_ifs.tests
    # run_test test_varenv17.sub
    run_test test_exp3.sub
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
    # run_test test_histexp7.sub - fails in CI
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
    # run_test test_unicode3.sub
    # run_test test_set-x.tests - using sets doesn't work, and this calls on several other tests here
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
    run_test test_exp10.sub
    # run_test test_nquote5.tests
    # run_test test_new-exp10.sub
    # run_test test_array5.sub
    # run_test test_test1.sub
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
    # run_test test_set-e.tests
    # run_test test_histexp3.sub - history expansion
    # run_test test_varenv22.sub - error script
    run_test test_procsub1.sub
    run_test test_comsub3.sub
    # run_test test_history6.sub
    # run_test test_builtins.tests
    # run_test test_extglob1.sub - extended globbing shouldn't work
    # run_test test_func4.sub
    # run_test test_alias.tests
    # run_test test_set-e3.sub - calls set-e3a.sub
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
    run_test test_dollar-at-star2.sub
    run_test test_exec12.sub
    run_test test_source1.sub
    run_test test_array10.sub
    run_test test_redir9.sub
    # run_test test_quotearray2.sub
    # run_test test_redir8.sub
    # run_test test_quotearray3.sub
    # run_test test_array11.sub
    # run_test test_dollar-at-star3.sub
    # run_test test_exec13.sub
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
    # run_test test_posixexp7.sub
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
    # run_test test_exec11.sub - traps could catch other system stuff not expected to work
    run_test test_dollar-at-star1.sub
    run_test test_new-exp5.sub
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
    # run_test test_exec3.sub - trap can have differing behavior
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
    run_test test_ifs1.su
    run_test test_getopts7.sub
    # run_test test_coproc.tests - error script
    run_test test_source4.sub
    run_test test_new-exp3.sub
    run_test test_dollar-at-star7.sub
    # run_test test_array29.sub
    # run_test test_array15.sub
    # run_test test_exportfunc.tests
    # run_test test_array14.sub
    run_test test_dollar-at-star6.sub
    # run_test test_new-exp2.sub
    run_test test_mapfile1.sub
    run_test test_source5.sub
    run_test test_getopts6.sub
    run_test test_errors4.sub
    run_test test_printf1.sub
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
