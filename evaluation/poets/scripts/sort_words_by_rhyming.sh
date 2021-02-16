# Sort words in Genesis by rhyming order.
INPUT=${INPUT:-inputs/genesis}
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort | uniq -c | rev | sort | rev 
