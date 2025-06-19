#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

mkdir -p "$OUT"


for input in $(cat "$PASH_TOP/evaluation/benchmarks/nlp/1000-books.txt" | head -n ${ENTRIES})
do
    cat $IN$input | grep -c 'light.\*light'                                 > ${OUT}${input}.stdout0
    cat $IN$input | grep -c 'light.\*light.\*light'                         > ${OUT}${input}.stdout1
    cat $IN$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > ${OUT}${input}.stdout2
done

echo 'done';
# rm -rf ${OUT}
