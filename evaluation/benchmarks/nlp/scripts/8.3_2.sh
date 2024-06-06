#!/bin/bash 
# tag: find_anagrams.sh
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.3_2/}
ENTRIES=${ENTRIES:-1000}
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
for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}/${input}.out
done

echo 'done';
# rm -rf "$OUT"
