#!/bin/bash

BENCHMARK_DIR=$(cd "$(dirname "$0")" && pwd)
scripts_dir="${BENCHMARK_DIR}/scripts"

export BENCHMARK_CATEGORY="web-search"
export LC_ALL=C

suffix=""
for arg in "$@"; do
    case "$arg" in
        --small) suffix="_small" ;;
        --min)   suffix="_min"   ;;
    esac
done

export BENCHMARK_DIR
export IN="${BENCHMARK_DIR}/inputs/index${suffix}.txt"
export WIKI="${BENCHMARK_DIR}/inputs/articles${suffix}/"
export WEB_INDEX_DIR="${BENCHMARK_DIR}/input"

export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/web-search.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$IN")"

bash "$scripts_dir/web-search.sh"
echo $?
