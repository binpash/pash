# 1-syllable words 
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[ 12*]' < ${INPUT} | grep -i '^[^aeiou]*[aeiou][^aeiou]*$' | sort | uniq -c | sed 5q
