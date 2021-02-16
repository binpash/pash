# Count words given an input file
# $INPUT is the input file
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' <  ${INPUT} | sort | uniq -c
