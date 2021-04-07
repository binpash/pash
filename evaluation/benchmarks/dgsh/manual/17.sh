#!/bin/bash
# from: https://github.com/dspinellis/dgsh/blob/master/example/reorder-columns.sh
# tag: reorder columns
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}

rm -f a b c d
mkfifo a b c d
cat a  | cut -d , -f 5-6 - | cat > c & 
cat b  | cut -d , -f 2-4 - | cat > d &
cat ${IN}/gsquarterly_december-2020.csv| tee a b > /dev/null  &
paste -d , c d
rm a b c d
