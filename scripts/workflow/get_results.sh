export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
rm -rf log_results
mkdir log_results

stats() (
    test_results_dir=$1  
    grep "are identical" "$test_results_dir"/result_status |
        sed "s,^$PASH_TOP/,," > log_results/$2_passed.log
    cat log_results/$2_passed.log >> log_results/passed.log
    grep  "are not identical" "$test_results_dir"/result_status |
        sed "s,^$PASH_TOP/,," > log_results/$2_failed.log
    # if the file has data, append it
    if [ -s log_results/$2_failed.log ] 
    then
        cat log_results/$2_failed.log >> log_results/failed.log
    else
        # remove since it's empty
        rm log_results/$2_failed.log
    fi
    TOTAL_TESTS=$(cat "$test_results_dir"/result_status | wc -l)
    PASSED_TESTS=$(grep  "are identical" "$test_results_dir"/result_status | wc -l)
    echo "$2: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
)

echo "Below follow the identical outputs:"     > log_results/passed.log
echo "Below follow the non-identical outputs:" > log_results/failed.log
#
## intro tests
stats "$PASH_TOP/evaluation/intro/output" intro
#
## Interface Tests
stats "$PASH_TOP/evaluation/tests/interface_tests/output" interface
#
## compiler Tests
stats "${PASH_TOP}/evaluation/tests/results" compiler
