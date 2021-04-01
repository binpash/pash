# words with no vowels
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c
