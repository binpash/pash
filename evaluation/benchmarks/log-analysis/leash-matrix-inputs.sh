#!/bin/bash

# Input matrix for evaluation/benchmarks/nlp.
# Mirrors the original nlp run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export LEASH_ENTRIES=84
    SCRIPT_INPUT_WIDTH=(
        "pcaps_heavy.sh:pcap_data_heavy/:1"
        "nginx_heavy.sh:log_data_heavy/:1"
    )
}
