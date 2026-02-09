#!/bin/bash

# Input matrix for evaluation/benchmarks/covid-mts.
# Mirrors the original covid-mts run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "1.sh:in_tiny.csv:2"
        # "2.sh:in.csv:16"
        # "3.sh:in.csv:16"
        # "4.sh:in.csv:16"
    )
}
