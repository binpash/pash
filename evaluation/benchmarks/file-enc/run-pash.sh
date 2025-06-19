#!/bin/bash

time ENTRIES=16 IN=log-analysis/inputs/pcap_data_heavy/ OUT=file-enc/outputs/ SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w1 --parallel_pipelines --parallel_pipelines_limit 8 scripts/compress_files_heavy.sh