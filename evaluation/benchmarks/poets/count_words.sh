# Count words given an input file
# $INPUT is the input file
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c
