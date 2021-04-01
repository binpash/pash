# Bigrams that appear twice
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | awk '$1 == 2 {print $2, $3}' 
