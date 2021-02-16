INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/}
mkdir -p $PASH_TOP/evaluation/scripts/input/output
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/scripts/input/output}
mogrify -format jpg -path ${OUTPUT} -thumbnail 100x100 ${INPUT}*.jpg
