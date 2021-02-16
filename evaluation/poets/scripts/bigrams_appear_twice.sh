# Bigrams that appear twice
INPUT=${INPUT:-inputs/genesis.bigrams}
awk '$1 == 2 {print $2, $3}' ${INPUT}
