#!/bin/bash
# compress all files in a directory
IN=${IN:-/dependency_untangling/pcap_data/}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/compress}

mkdir -p ${OUT}

for item in $(hdfs dfs -ls -C ${IN});
do
    output_name=$(basename $item).zip
    hdfs dfs -cat $item | gzip -c > $OUT/$output_name
done

echo 'done';
