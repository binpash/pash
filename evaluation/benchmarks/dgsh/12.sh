#!/bin/bash
#
# SYNOPSIS Parallel word count
# DESCRIPTION
# Number of processes
# word count
#FIXME dataset
set -e
N=8
# Collation order for sorting
export LC_ALL=C
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input/dblp.xml}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}

# split our chunk in N files
split -dn ${IN}/$N $1
# our files are in x00 x01 x02 
for i in $(seq $((N -1))) ; 
do 
    f=x0$i
    cat $f | tr -s ' \t\n\r\f' '\n' | sort -S 512M  > ${OUT}/${f}.res
done 

# i don't like it at all, u def can optimize it
cat ${OUT}/*.res | sort -S 512M  | uniq -c
