#!/bin/bash
# encrypt all files in a directory 
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/encrypt}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
run_tests() {
    openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -in $1 -out $OUT/$(basename $1).enc -k 'key'
}

export -f run_tests
pkg_count=0

for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > ${LOGS}/${pkg_count}.log
done

echo 'done';
