# Merge upper and lower counts
# ${INPUT} is the input file
# ${OUTPUT} is the generated file
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
tr '[a-z]' '[A-Z]' < ${INPUT} | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c
