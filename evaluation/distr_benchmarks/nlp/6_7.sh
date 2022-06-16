#!/bin/bash
# tag: verse_2om_3om_2instances
# set -e
# verses with 2 or more, 3 or more, exactly 2 instances of light.

IN=${IN:-/nlp/pg/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/nlp/output/6_7/}
ENTRIES=${ENTRIES:-1060}
mkdir -p "$OUT"

for input in $(hdfs dfs -ls -C ${IN} | head -n ${ENTRIES} | xargs -n 1 -I arg1 basename arg1)
do
    hdfs dfs -cat -ignoreCrc $IN/$input | grep -c 'light.\*light'                                 > ${OUT}/${input}.out0
    hdfs dfs -cat -ignoreCrc $IN/$input | grep -c 'light.\*light.\*light'                         > ${OUT}/${input}.out1
    hdfs dfs -cat -ignoreCrc $IN/$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > ${OUT}/${input}.out2
done

echo 'done';
rm -rf ${OUT}
