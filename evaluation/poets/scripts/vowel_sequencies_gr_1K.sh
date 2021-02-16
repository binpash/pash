# Vowel sequencies that appear >=1K times
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk '$1 >= 1000'
