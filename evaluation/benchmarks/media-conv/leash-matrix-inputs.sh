#!/bin/bash

# Input matrix for evaluation/benchmarks/nlp.
# Mirrors the original nlp run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export ENTRIES=150
    SCRIPT_INPUT_WIDTH=(
        "to_mp3_heavy.sh:wav_full_heavy/wav/:1",
    )
}
