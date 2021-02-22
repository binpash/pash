INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/posh/log20170314.csv}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/scripts/input/posh/log20170315.csv}
cat ${INPUT} ${INPUT2} | grep 'bar' > $PASH_TOP/evaluation/scripts/input/output.txt
