#!/bin/bash

# Input matrix for evaluation/benchmarks/max-temp.
# Mirrors the original max-temp run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "max-temp-process.sh:temperatures.2015.txt:16"
        "max-temp-process.sh:temperatures.2015.txt:32"
        "max-temp-process.sh:temperatures.2015.txt:64"
        # "max-temp-process.sh:temperatures_1G.txt:16"
    )
}
