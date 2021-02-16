# count consonant sequences
INPUT=${INPUT:-inputs/genesis}
tr '[a-z]' '[A-Z]' < ${INPUT} | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c
