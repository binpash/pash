#!/bin/bash

# Input matrix for evaluation/benchmarks/oneliners.
# Mirrors the original oneliners input selections (including commented candidates).
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        # "sort.sh:1G.txt:4"
        # "sort-sort.sh:3G.txt:16"
        "sort.sh:1M.txt:2"
        # "sort-sort.sh:3G.txt:16"
        # "wf.sh:3G.txt:16"
        # "top-n.sh:3G.txt:16"
        # "set-diff.sh:3G.txt:16"
        # "bi-grams.sh:1G.txt:16"
        # "spell.sh:3G.txt:16"
        # "nfa-regex.sh:3G.txt:128"
    )
}
