#!/bin/bash
# tag: trigram_rec
set -e

# FIXME: what is this?
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/input/output/}
TER=${TER:-$PASH_TOP/pa.sh}


trigrams() {
  ls $IN/ | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.words
  tail +2 ${OUT}.words > ${OUT}.nextwords
  tail +3 ${OUT}.words > ${OUT}.nextwords2
  paste ${OUT}.words ${OUT}.nextwords ${OUT}.nextwords2 |
  sort | uniq -c  > ${OUT}.trigrams
}

# Recursive PaSh calls to trigrams (defined below)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
# Store it because we call other script with different INPUT


# we should pipe the ouput of grep to the input of pash. Since pash does not
# offer a way to accept command line args to alter the hardcoded input
# we set the env variable as the commented out code below
# Old 
#grep 'the land of' $IN | $PASH_TOP/pa.sh scripts/count_trigrams.sh | sort -nr | sed 5q
#grep 'And he said' $IN | $PASH_TOP/pa.sh scripts/count_trigrams.sh | sort -nr | sed 5q


ls ${IN} | sed "s;^;$IN;"| xargs cat | grep 'the land of' | ${TER} 4_3b.sh | sort -nr | sed 5q
ls ${IN} | sed "s;^;$IN;"| xargs cat | grep 'And he said' | ${TER} 4_3b.sh | sort -nr | sed 5q



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
