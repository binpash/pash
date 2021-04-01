# Uppercase words by token
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^[A-Z]'
#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^[A-Z]'| paste -sd+ | bc
