#!/bin/bash
# tag: generate thumbnails with size 100x100 
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/images/originals}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/posh/input/output}
mogrify  -format gif -path ${OUT} -thumbnail 100x100 ${IN}/*.jpg
