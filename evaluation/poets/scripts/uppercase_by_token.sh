# Uppercase words by token
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^[A-Z]'| paste -sd+ | bc
