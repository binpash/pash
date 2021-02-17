# Recursive PaSh calls to trigrams (defined below)
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
# Store it because we call other script with different INPUT
IN=${INPUT}
GREP_RES=${IN}.rec
grep 'the land of' ${INPUT} > ${GREP_RES}
export INPUT=${GREP_RES}
$PASH_TOP/pa.sh scripts/count_trigrams.sh | cat ${GREP_RES}.trigrams | sort -nr | sed 5q
# reset the input to genesis
export INPUT=${IN}
# process on genesis
grep 'And he said' ${INPUT} > ${GREP_RES}2
export INPUT=${GREP_RES}2
bash scripts/count_trigrams.sh | cat ${GREP_RES}2.trigrams | sort -nr | sed 5q
