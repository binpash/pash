#!/bin/bash
# tag: find_anagrams.sh
# set -e

mkdir -p "$OUT"


pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    sort -u > ${TEMPDIR}/${input}.types
    rev < ${TEMPDIR}/${input}.types > ${TEMPDIR}/${input}.types.rev
    sort ${TEMPDIR}/${input}.types ${TEMPDIR}/${input}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}"
    rm -rf ${TEMPDIR}
}

export -f pure_func
for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
