# Find anagrams 
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
# need to generate words
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} > ${INPUT}.words
# need to generate types
sort -u ${INPUT}.words > ${INPUT}.types
# Actual find anagram script
rev < ${INPUT}.types > ${INPUT}.types.rev
sort ${INPUT}.types ${INPUT}.types.rev |
uniq -c |
awk '$1 >= 2 {print $2}'
