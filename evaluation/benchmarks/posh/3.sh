#!/bin/bash
# tag: ray tracing
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/ray_tracing/logs}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/posh/input/output}
echo $IN
cat ${IN}/1.INFO | grep "\[RAY\]" | head -n1 | cut -c 7- > ${OUT}/rays.csv
cat ${IN}/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- >> ${OUT}/rays.csv
cat ${OUT}/rays.csv | sed -n '/^590432,/p' > ${OUT}/output.log
