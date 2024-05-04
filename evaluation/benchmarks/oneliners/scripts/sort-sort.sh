#!/bin/bash
# Calculate sort twice

cd "$(dirname "$0")" || exit 1

SIZE=500M

IN="input/$SIZE.txt"

cat $IN | tr A-Z a-z | sort | sort -r
