#!/bin/bash 
# tag: bigrams_appear_twice.sh
# set -e

# Calculate the bigrams (based on 4_3.sh script)
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
#ls ${IN} | sed "s;^;$IN;"| xargs cat | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}.input.words
#tail +2 ${OUT}.input.words > ${OUT}.input.nextwords
#paste ${OUT}.input.words ${OUT}.input.nextwords | sort | uniq -c > ${OUT}.input.bigrams
## find the bigrams that appear exactly twice 
#awk "\$1 == 2 {print \$2, \$3}" ${OUT}.input.bigrams
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/8.2_2}
mkdir -p "$OUT"
ENTRIES=${ENTRIES:-1000}
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | tr -sc '[A-Z][a-z]' '[\012*]' > "${OUT}/$(basename ${input}).input.words"
    tail +2 "${OUT}/$(basename ${input}).input.words" > "${OUT}/$(basename ${input}).input.nextwords"
    paste "${OUT}/$(basename ${input}).input.words" "${OUT}/$(basename ${input}).input.nextwords" | sort | uniq -c > "${OUT}/$(basename ${input}).input.bigrams"
    awk "\$1 == 2 {print \$2, \$3}" "${OUT}/$(basename ${input}).input.bigrams"
done

for output in $(ls ${OUT} | sed "s;^;$OUT/;")
do
    cat $output
done

echo 'done';
rm -rf "${OUT}"
