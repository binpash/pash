# words with no vowels
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -vi '[aeiou]' | sort | uniq -c
