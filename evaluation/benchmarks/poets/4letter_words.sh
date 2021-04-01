# four-letter words
# the original script has both versions
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' 
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$' 

#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | grep -c '^....$' | paste -sd+ | bc
#tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT} | sort -u | grep -c '^....$' | paste -sd+ | bc

