#!/bin/bash
# tag: bigrams.sh
# set -e

# Bigrams (contrary to our version, this uses intermediary files)
mkdir -p "$OUT"


pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.input.words
    tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
    paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords
    rm -rf ${TEMPDIR}
}
export -f pure_func

for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$"| pure_func $input| sort | uniq -c > ${OUT}${input}.input.bigrams
done

echo 'done';
# rm -rf ${OUT}
