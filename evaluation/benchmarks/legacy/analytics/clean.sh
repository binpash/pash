#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

TOP=$(git rev-parse --show-toplevel)
input_dir="${TOP}/evaluation/benchmarks/analytics/inputs"
outputs_dir="${TOP}/evaluation/benchmarks/analytics/outputs"

rm -rf "$outputs_dir"

if [ "$force" = true ]; then
    rm -rf "$input_dir"
    rm -rf "${TOP}/evaluation/benchmarks/analytics/go_install"
fi
