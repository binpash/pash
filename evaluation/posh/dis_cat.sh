INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/log1}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/scripts/input/log2}
cat ${INPUT} ${INPUT2} | grep 'bar' > $PASH_TOP/evaluation/scripts/input/output.txt
