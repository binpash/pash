# Sort words in Genesis by rhyming order.
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev 
