# Count tri-grams
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} > ${INPUT}.words
tail +2 ${INPUT}.words > ${INPUT}.nextwords
tail +3 ${INPUT}.words > ${INPUT}.nextwords2
paste ${INPUT}.words ${INPUT}.nextwords ${INPUT}.nextwords2 |
sort | uniq -c > ${INPUT}.trigrams
