#!/bin/bash
# tag: bigrams_appear_twice.sh
# set -e

# Calculate the bigrams (based on 4_3.sh script)
mkdir -p "$OUT"


pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.input.words
    tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
    paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords | sort | uniq -c > ${TEMPDIR}/${input}.input.bigrams
    awk "\$1 == 2 {print \$2, \$3}" ${TEMPDIR}/${input}.input.bigrams
    rm -rf {TEMPDIR}
}

export -f pure_func
for input in $(echo $inputs | tr " " "\n" | head -n ${ENTRIES})
do
    cat $IN$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > ${OUT}${input}.stdout
done

echo 'done';
# rm -rf "$OUT"
