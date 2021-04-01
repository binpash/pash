# Bigrams that appear twice
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis.bigrams}
ls $IN/ | xargs cat | awk '$1 == 2 {print $2, $3}' 
