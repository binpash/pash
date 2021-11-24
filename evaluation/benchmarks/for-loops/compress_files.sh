#!/bin/bash
# compress all files in a directory
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/pcaps/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/output/compress}
LOGS=${OUT}/logs
mkdir -p ${OUT}/logs
run_tests() {
    echo ${OUT}
    name=$(basename $1).zip
    zip -r ${OUT}/$name $1
}

export -f run_tests

pkg_count=0
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > ${LOGS}/$pkg_count.log
done
