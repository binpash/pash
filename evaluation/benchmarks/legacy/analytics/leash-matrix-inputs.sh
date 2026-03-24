#!/bin/bash

# Input matrix for evaluation/benchmarks/nlp.
# Mirrors the original nlp run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export LEASH_ENTRIES=84
    SCRIPT_INPUT_WIDTH=(
        # "nginx_heavy.sh:nginx-logs_full/:1" #TODO FIX
        # "pcaps_heavy.sh:pcaps_full/:1"
        # "nginx.sh:nginx-logs_small:1"
        # "pcaps.sh:pcaps_small:1"
        "port-scan.sh:port_scan_small:1"
        # "ray-tracing.sh:ray_tracing_small:1"
    )
}
