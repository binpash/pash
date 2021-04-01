# Merge upper and lower counts
# ${INPUT} is the input file
# ${OUTPUT} is the generated file
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr '[a-z]' '[A-Z]' | tr -sc '[A-Z]' '[\012*]' | sort | uniq -c
