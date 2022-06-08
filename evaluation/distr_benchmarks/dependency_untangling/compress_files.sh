#!/bin/bash
# compress all files in a directory
IN=${IN:-/dependency_untangling/pcap_data/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/compress}

mkdir -p ${OUT}
pure_func() {
    zip -r --
}

export -f pure_func

for item in $(hdfs dfs -ls -C ${IN});
do
    output_name=$(basename $item).zip
    hdfs dfs -cat $item | pure_func > $OUT/$output_name
done

echo 'done';
