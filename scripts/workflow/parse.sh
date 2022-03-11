bash ./get_results.sh > out
mv out log_results
cat log_results/out
while read p; do
    PASSED=$(echo $p |  awk -F'[^0-9]+' '{ print $2 }')
    TOTAL=$(echo $p |  awk -F'[^0-9]+' '{ print $3 }')
    FAILED=$((passed - failed))
    # failed, print to stdout
    if [ $PASSED -ne $TOTAL ]; then
      # get the benchmark name
      f=${p%% *}
      # strip the :
      f="${f%?}"
      # dump the failed tests
      cat log_results/${f}_failed.log
    fi
done < log_results/out
