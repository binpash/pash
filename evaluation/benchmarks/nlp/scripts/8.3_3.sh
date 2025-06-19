#!/bin/bash
# tag: compare_exodus_genesis.sh
# set -e

mkdir -p "$OUT"


pure_func() {
    input=$1
    input2=$2
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}1.types
    cat  ${input2} | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${TEMPDIR}/${input}2.types
    sort ${TEMPDIR}/${input}1.types ${TEMPDIR}/${input}2.types ${TEMPDIR}/${input}2.types | uniq -c | head
    rm -rf ${TEMPDIR}
}
export -f pure_func

for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | pure_func $input $INPUT2 > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
