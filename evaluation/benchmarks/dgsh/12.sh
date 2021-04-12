#!/bin/bash
# from https://www2.dmst.aueb.gr/dds/sw/dgsh/#parallel-word-count
# tag: parallel word count
#FIXME dataset, outdated split ?
set -e
N=8
# Collation order for sorting
export LC_ALL=C
IN=${FULL:-$PASH_TOP/evaluation/benchmarks/dgsh/input/dblp.xml}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}

# split our chunk in N files
split -dn  $N ${IN}
# our files are in x00 x01 x02 
for i in $(seq $((N -1))) ; 
do 
    f=x0$i
    cat $f | tr -s ' \t\n\r\f' '\n' | sort -S 512M  > ${OUT}/${f}.res
done 

# i don't like it at all, u def can optimize it
cat ${OUT}/*.res | sort -S 512M  | uniq -c
