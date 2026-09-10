#!/bin/bash

# Input matrix for evaluation/benchmarks/agent (log-severity summary).
# The shared runner uses IN="agent/inputs/<INPUT>" and OUT on S3, so INPUT is the
# S3 key prefix of the log directory: "logs/".  example.sh iterates logs.txt.
set_leash_benchmark_inputs() {
    SCRIPT_INPUT_WIDTH=(
        "terminal-bench-log-summary.sh:logs/:1"
        # "example.sh:logs/:2"
        # "example.sh:logs/:4"
    )
}
