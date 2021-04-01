# Sort words in Genesis by rhyming order.
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev 
