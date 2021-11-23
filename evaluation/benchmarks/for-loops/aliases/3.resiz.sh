#!/bin/bash

# tag: resize image 
set -e
IN=${JPG:-$PASH_TOP/evaluation/benchmarks/aliases/input/jpg}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/aliases/input/out}
# could we convert it to somelthing like this
rm -f $OUT/*
for i in $IN/*.jpg; 
do 
out=$OUT/$(basename -- $i)
convert -resize 70% "$i" "$out"; 
done
# or this is good for us
#find $IN -name "*.jpg" | 
#  xargs -n1 basename |
#  sed "s;\(.*\);-resize 70% $IN/\1 $OUT/\1.70;" |
#  xargs -L1  convert
