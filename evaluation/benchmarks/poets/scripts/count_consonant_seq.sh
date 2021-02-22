# count consonant sequences
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
tr '[a-z]' '[A-Z]' < ${INPUT} | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c
