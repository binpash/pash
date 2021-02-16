# Count vowel Sequences
INPUT=${INPUT:-inputs/genesis}
tr 'a-z' '[A-Z]' < ${INPUT} | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c 
