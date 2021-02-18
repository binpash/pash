# Uppercase words by type
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^[A-Z]'
#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^[A-Z]'|  paste -sd+ | bc
