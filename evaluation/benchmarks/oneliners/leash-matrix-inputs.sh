#!/bin/bash

# Input matrix for evaluation/benchmarks/oneliners.
# Mirrors the original oneliners input selections (including commented candidates).
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "sort.sh:3G.txt:64"
        "sort-sort.sh:3G.txt:64"
        "wf.sh:3G.txt:64"
        "top-n.sh:3G.txt:64"
        "set-diff-leash.sh:3G.txt:64"
        "bi-grams.sh:3G.txt:64"
        "spell.sh:3G.txt:64"
        "nfa-regex.sh:3G.txt:64"
    )
}
