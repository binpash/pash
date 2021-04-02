#!/bin/bash
# works on pash
# Adapted from the DGSH
# https://github.com/dspinellis/dgsh/blob/master/example/reorder-columns.sh
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input/gsquarterly_december-2020.csv}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}

set -e
#OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/dgsh/output}
#cd ${OUTPUT}
rm -f a b c d
mkfifo a b c d
cat a  | cut -d , -f 5-6 - | cat > c & 
cat b  | cut -d , -f 2-4 - | cat > d &
cat ${IN} | tee a b > /dev/null  &
paste -d , c d
rm a b c d
