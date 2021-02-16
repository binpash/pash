# words with no vowels
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -vi '[aeiou]' | sort | uniq -c
