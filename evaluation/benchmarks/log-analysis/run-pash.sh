#!/bin/bash

time ENTRIES=16 IN=log-analysis/inputs/pcap_data_heavy/ OUT=log-analysis/outputs/pcap-heavy-pash-core4/ SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w1 --parallel_pipelines --parallel_pipelines_limit 8 scripts/pcaps_heavy.sh

# time ENTRIES=84 IN=log-analysis/inputs/log_data_heavy/ OUT=log-analysis/outputs/nginx-heavy-pash-core4/ SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w1 --parallel_pipelines --parallel_pipelines_limit 4 scripts/nginx_heavy.sh
