#!/bin/bash
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/node_modules}
MIR_BIN=${MIR_BIN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/mir-sa/.bin/mir-sa}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/mir}
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

echo 'done';
