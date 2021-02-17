# Sort words in Genesis by folding case.
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort | uniq -c | sort -f 
