#!/bin/bash
# tag: wav-to-mp3
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/wav}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/mp3}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
run_tests(){
    FILE=$1
    ffmpeg -y -i $FILE -f mp3 -ab 192000 $OUT/$(basename $FILE).mp3  2>/dev/null
}

export -f run_tests

pkg_count=0
for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > ${LOGS}/${pkg_count}.log
done

echo 'done';
