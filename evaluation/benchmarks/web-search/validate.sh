#!/bin/bash

BENCHMARK_DIR=$(cd "$(dirname "$0")" && pwd)

status=0
for f in 1-grams.txt 2-grams.txt; do
    if [ ! -s "$BENCHMARK_DIR/$f" ]; then
        echo "FAIL: $f missing or empty" >&2
        status=1
    fi
done

if [ $status -eq 0 ]; then
    echo "web-search 0"
fi
exit $status
