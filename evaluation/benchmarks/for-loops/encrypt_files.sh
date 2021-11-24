#!/bin/bash
# encrypt all files in a directory 
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/encrypt}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
run_tests() {
    echo 'key' | openssl enc -aes-256-cbc -e -md sha512 -in $item -out $OUT/$(basename $1).enc --pass stdin 
}

export -f run_tests
pkg_count=0

for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item  > ${LOGS}/${pkg_count}.log
done
