#!/bin/bash

# tag: resize image 
set -e

IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/jpg
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
# do we need such a catch-all find if we know there's only jpg?
# find . -name "*.jpg" -o -name "*.png" -o -name "*.jpeg" -o -name "*.JPG" -o -name "*.PNG" -o -name "*.JPEG" 
find $IN -name "*.jpg" -o -name "*.jpeg" -o -name "*.JPG" -o  -name "*.JPEG" | 
  xargs -n1 basename |
  sed "s;\(.*\);$IN/\1 $OUT/\1.70;" |
  tr '\n' '\0' |
  xargs -0 -n1 echo convert -resize 70% 
