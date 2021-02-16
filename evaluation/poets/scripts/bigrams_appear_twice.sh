# Bigrams that appear twice
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis.bigrams}
awk '$1 == 2 {print $2, $3}' ${INPUT}
