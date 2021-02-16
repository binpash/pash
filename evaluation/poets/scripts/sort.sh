# Sort
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort | uniq -c | sort -nr

