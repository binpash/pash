#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_7/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    cat $IN/$input | grep -c 'light.\*light'                                 > ${OUT}/${input}.out0
    cat $IN/$input | grep -c 'light.\*light.\*light'                         > ${OUT}/${input}.out1
    cat $IN/$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > ${OUT}/${input}.out2
done

echo 'done';
# rm -rf ${OUT}
