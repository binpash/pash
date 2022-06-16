#!/bin/bash
# encrypt all files in a directory 
IN=${IN:-/dependency_untangling/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/encrypt}
mkdir -p ${OUT}

pure_func() {
    openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k 'key'
}
export -f pure_func

for item in $(hdfs dfs -ls -C ${IN});
do
    output_name=$(basename $item).enc
    hdfs dfs -cat -ignoreCrc $item | pure_func > $OUT/$output_name
done

echo 'done';
