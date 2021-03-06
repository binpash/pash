# Count vowel Sequences
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
tr 'a-z' '[A-Z]' < ${INPUT} | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c 
