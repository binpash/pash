# four-letter words
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^....$' | paste -sd+ | bc
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^....$' | paste -sd+ | bc

