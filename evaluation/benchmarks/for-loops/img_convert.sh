#!/bin/bash
# tag: resize image 
set -e
IN=${JPG:-$PASH_TOP/evaluation/benchmarks/for-loops/input/jpg}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/jpg}
mkdir -p ${OUT}
for i in $IN/*.jpg; 
do 
    out=$OUT/$(basename -- $i)
    convert -resize 70% "$i" "$out"; 
done
