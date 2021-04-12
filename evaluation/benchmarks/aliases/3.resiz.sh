#!/bin/bash

# tag: resize image 
set -e


IN=${JPG:-$PASH_TOP/evaluation/benchmarks/aliases/input/jpg}
OUT=${OUT:-PASH_TOP/evaluation/benchmarks/aliases/input/out}

find $IN -name "*.jpg" | 
  xargs -n1 basename |
  sed "s;\(.*\);-resize 70% $IN/\1 $OUT/\1.70;" |
  xargs -L1  convert
