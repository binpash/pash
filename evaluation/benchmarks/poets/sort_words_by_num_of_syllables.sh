# Sort words by number of syllables
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${INPUT}.words
tr -sc '[AEIOUaeiou\012]' ' ' < ${INPUT}.words | awk '{print NF}' > ${INPUT}.syl
paste ${INPUT}.syl ${INPUT}.words | sort -nr | sed 5q
