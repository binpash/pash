#!/bin/bash
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/node_modules}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/output/mir}
MIR_BIN=${IN}/.bin/mir-sa
mkdir -p ${OUT}/
pkg_count=0
run_tests() {
    cd $1;
    ${MIR_BIN} -p   2>>${OUT}/error.log
}
export -f run_tests
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > ${OUT}/$pkg_count.log
done
