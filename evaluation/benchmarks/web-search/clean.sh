#!/bin/bash

BENCHMARK_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$BENCHMARK_DIR" || exit 1

rm -f 1-grams.txt 2-grams.txt 3-grams.txt

# Remove any leftover named pipes
find . -maxdepth 1 -type p -delete
