# Sort words by number of syllables
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u > ${INPUT}.words
tr -sc '[AEIOUaeiou\012]' ' ' < ${INPUT}.words | awk '{print NF}' > ${INPUT}.syl
paste ${INPUT}.syl ${INPUT}.words | sort -nr | sed 5q
