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
for input in $(ls ${IN} | head -n ${ENTRIES} | sed "s;^;$IN;")
do
    cat "$input" | grep -c 'light.\*light'                                 > "${OUT}/$(basename ${input})"
    cat "$input" | grep -c 'light.\*light.\*light'                         > "${OUT}/$(basename ${input})1"
    cat "$input" | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > "${OUT}/$(basename ${input})2"
done

for output in $(ls ${OUT} | sed "s;^;$OUT;")
do
    cat $output
done

echo 'done';
rm -rf ${OUT}
