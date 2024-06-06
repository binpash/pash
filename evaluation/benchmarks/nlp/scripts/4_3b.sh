#!/bin/bash 
#tag: count_trigrams.sh
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3b/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.words
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords2
    paste ${TEMPDIR}/${input}.words ${TEMPDIR}/${input}.nextwords ${TEMPDIR}/${input}.nextwords2 |
    sort | uniq -c 
    rm -rf ${TEMPDIR}
}
export -f pure_func
for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}/${input}.trigrams
done

echo 'done';
# rm -rf ${OUT}
