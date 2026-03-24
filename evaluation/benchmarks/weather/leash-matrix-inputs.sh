#!/bin/bash

# Input matrix for evaluation/benchmarks/weather.
# Mirrors the original weather run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "temp-analytics.sh:temperatures.2015.txt:64"
    )
}
