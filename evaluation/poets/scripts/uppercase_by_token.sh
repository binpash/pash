# Uppercase words by token
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^[A-Z]'
