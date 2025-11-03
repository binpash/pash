#!/bin/bash

cd $(dirname $0)

mkdir -p outputs

{ time ENTRIES=1 IN=log-analysis/inputs/pcap_data_heavy/ OUT=file-enc/outputs/ SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w1 --parallel_pipelines --parallel_pipelines_limit 4 scripts/compress_files_heavy.sh; } &>outputs/localecorrect-compress_files_heavy-pash-core4-time.log

sleep 20
logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
if [ -d "$logs_dir" ]; then
    echo "Removing existing logs directory: $logs_dir"
    rm -rf "$logs_dir"
fi
python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"