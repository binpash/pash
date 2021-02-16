# Uppercase words by token
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^[A-Z]'
