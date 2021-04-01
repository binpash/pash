#!/bin/bash

# tag: resize image 
set -e

IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/jpg
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
#find $IN -name "*.jpg" | 
#  xargs -n1 basename |
#  sed "s;\(.*\);$IN/\1 $OUT/\1.70;" |
#  tr '\n' '\0' |
#  xargs -0 -n1 convert -resize 70% 
# can be optimized 
find . -iname "*.jpg" -printf "%p  ${OUT}/%f.70\n" |
    xargs -r -n2 convert -resize 70%
