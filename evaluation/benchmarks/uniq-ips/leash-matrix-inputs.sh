#!/bin/bash

# Input matrix for evaluation/benchmarks/uniq-ips.
# Mirrors the original uniq-ips input selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "uniq-ips.sh:logs-popcount-org.txt:16"
        "uniq-ips.sh:logs-popcount-org.txt:32"
        "uniq-ips.sh:logs-popcount-org.txt:64"
        # Keep candidate inputs here if more cases are added later.
    )
}
