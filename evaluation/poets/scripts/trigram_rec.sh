OUT=$PASH_TOP/evaluation/script/input
trigrams()
(
    tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} > ${INPUT}.words
    tail +2 ${INPUT}.words > ${INPUT}.nextwords
    tail +3 ${INPUT}.words > ${INPUT}.nextwords2
    paste ${INPUT}.words ${INPUT}.nextwords ${INPUT}.nextwords2 |
    sort | uniq -c  > ${INPUT}.trigrams
)

# Recursive PaSh calls to trigrams (defined below)
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
# Store it because we call other script with different INPUT


# we should pipe the ouput of grep to the input of pash. Since pash does not
# offer a way to accept command line args to alter the hardcoded input
# we set the env variable as the commented out code below
grep 'the land of' ${INPUT} | $PASH_TOP/pa.sh scripts/count_trigrams.sh | sort -nr | sed 5q
grep 'And he said' ${INPUT} | $PASH_TOP/pa.sh scripts/count_trigrams.sh | sort -nr | sed 5q



#IN=${INPUT}
#GREP_RES=${IN}.rec
#grep 'the land of' ${INPUT} > ${GREP_RES}
#export INPUT=${GREP_RES}
#$PASH_TOP/pa.sh scripts/count_trigrams.sh | cat ${GREP_RES}.trigrams | sort -nr | sed 5q
## reset the input to genesis
#export INPUT=${IN}
## process on genesis
#grep 'And he said' ${INPUT} > ${GREP_RES}2
#export INPUT=${GREP_RES}2
#$PASH_TOP/pa.sh scripts/count_trigrams.sh | cat ${GREP_RES}2.trigrams | sort -nr | sed 5q
