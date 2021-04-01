# Compare Exodus and Genesis
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/scripts/input/exodus}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u >  1.types
tr -sc '[A-Z][a-z]' '[\012*]' < ${INPUT2} | sort -u > 2.types
# do we really need the same thing twice
sort ${INPUT}.types ${INPUT2}.types ${INPUT2}.types | uniq -c | head
