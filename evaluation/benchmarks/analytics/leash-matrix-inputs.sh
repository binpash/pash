#!/bin/bash

# Input matrix for evaluation/benchmarks/nlp.
# Mirrors the original nlp run-leash.sh selection.
# Each row is script:input:width[:entries].
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export LEASH_ENTRIES=1
    SCRIPT_INPUT_WIDTH=(
        # "pcaps_heavy.sh:pcap_data_heavy/:1:16"
        # "nginx_heavy.sh:log_data_heavy/:1:84"
        # "encrypt_files.sh:pcap_data_heavy/:1:16"
        "img_convert.sh:jpg_full/jpg/:1:1"
    )
}
