#!/bin/bash
# tag: resize image 
IN=${JPG:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/jpg}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/jpg}
mkdir -p ${OUT}
for i in $IN/*.jpg; 
do 
    out=$OUT/$(basename -- $i)
    convert -resize 70% "$i" "$out"; 
done

echo 'done';
