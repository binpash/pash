#!/bin/bash
# tag: wav-to-mp3
IN=${IN:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/wav}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/mp3}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
pure_func(){
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1  2>/dev/null
}

export -f pure_func

pkg_count=0
for item in $(hdfs dfs -ls -C /for-loops/wav);
do
    pkg_count=$((pkg_count + 1));
    out="$OUT/$(basename $item).mp3"
    hdfs dfs -cat $item | pure_func > $out
done

echo 'done';
