# Uppercase words by type
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^[A-Z]'
