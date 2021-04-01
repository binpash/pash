# Find anagrams 
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
# need to generate words
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${INPUT}.words
# need to generate types
sort -u ${INPUT}.words > ${INPUT}.types
# Actual find anagram script
rev < ${INPUT}.types > ${INPUT}.types.rev
sort ${INPUT}.types ${INPUT}.types.rev |
uniq -c |
awk '$1 >= 2 {print $2}'
