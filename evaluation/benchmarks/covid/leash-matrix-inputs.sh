#!/bin/bash

# Input matrix for evaluation/benchmarks/covid.
# Mirrors the original covid run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "1.sh:in.csv:64"
        "2.sh:in.csv:64"
        "3.sh:in.csv:64"
        "4.sh:in.csv:64"
    )
}
