# Sort words in Genesis by rhyming order.
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort | uniq -c | rev | sort | rev 
