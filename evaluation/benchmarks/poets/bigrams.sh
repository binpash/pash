# Bigrams (contrary to our version, this uses intermediary files)
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > input.words
tail +2 input.words > input.nextwords
paste input.words input.nextwords | sort | uniq -c > input.bigrams
