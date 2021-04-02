# works on pash
# Adapted from the DGSH
# https://github.com/dspinellis/dgsh/blob/master/example/reorder-columns.sh
INPUT=${INPUT:-$PASH_TOP/evaluation/dgsh/input/gsquarterly_december-2020.csv}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/dgsh/output}
cd ${OUTPUT}

mkfifo a b c d
cat a  | cut -d , -f 5-6 - | cat > c & 
cat b  | cut -d , -f 2-4 - | cat > d &
cat ${INPUT} | tee a b > /dev/null  &
paste -d , c d
rm a b c d
