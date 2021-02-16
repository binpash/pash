# words with no vowels
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -vi '[aeiou]' | sort | uniq -c
