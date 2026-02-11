#!/bin/bash

# Input matrix for evaluation/benchmarks/nlp.
# Mirrors the original nlp run-leash.sh selection.
# Size flags are currently ignored for this benchmark.
set_leash_benchmark_inputs() {
    export ENTRIES=1000
    SCRIPT_INPUT_WIDTH=(
        # "1_1.sh:pg_heavy/:1"
        # "2_1.sh:pg_heavy/:1"
        # "2_2.sh:pg_heavy/:1"
        # "3_1.sh:pg_heavy/:1"
        # "3_2.sh:pg_heavy/:1"
        # "3_3.sh:pg_heavy/:1"
        # "4_3.sh:pg_heavy/:1"
        # "6_1_1.sh:pg_heavy/:1"
        # "6_1_2.sh:pg_heavy/:1"
        # "6_2.sh:pg_heavy/:1"
        # "6_3.sh:pg_heavy/:1"
        "6_4.sh:pg_heavy/:1"
        # "6_5.sh:pg_heavy/:1"
        # "6_7.sh:pg_heavy/:1"
        # "7_1.sh:pg_heavy/:1"
        # "7_2.sh:pg_heavy/:1"
        # "8.2_1.sh:pg_heavy/:1"
        # "6_1.sh:pg_heavy/:1"
        # "8_1.sh:pg_heavy/:1"
        # "8.3_2.sh:pg_heavy/:1"
        # "8.3_3.sh:pg_heavy/:1"
    )
}
