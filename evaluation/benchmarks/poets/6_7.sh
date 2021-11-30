#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/poets/output/6_7/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"
#ls ${IN} | sed "s;^;$IN;"| xargs cat | grep -c 'light.\*light'  
#ls ${IN} | sed "s;^;$IN;"| xargs cat | grep -c 'light.\*light.\*light' 
#ls ${IN} | sed "s;^;$IN;"| xargs cat | grep 'light.\*light' | grep -vc 'light.\*light.\*light'
for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat "$IN/$input" | grep -c 'light.\*light'                                 > "${OUT}/${input}.out0"
    cat "$IN/$input" | grep -c 'light.\*light.\*light'                         > "${OUT}/${input}.out1"
    cat "$IN/$input" | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > "${OUT}/${input}.out2"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf ${OUT}
