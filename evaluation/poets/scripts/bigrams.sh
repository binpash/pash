# Bigrams (contrary to our version, this uses intermediary files)
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} > ${INPUT}.words
tail +2 ${INPUT}.words > ${INPUT}.nextwords
paste ${INPUT}.words ${INPUT}.nextwords |
sort | uniq -c > ${INPUT}.bigrams
