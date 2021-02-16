# four-letter words
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^....$'
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^....$'
