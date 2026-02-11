#!/bin/bash

# Input matrix for evaluation/benchmarks/file-enc.
# Mirrors the original file-enc run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export ENTRIES=16
    SCRIPT_INPUT_WIDTH=(
        "compress_files_heavy.sh:pcap_data_heavy/:1"
    )
}
