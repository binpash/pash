#!/bin/bash
# Calculate sort twice

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

cat $IN | tr A-Z a-z | sort | sort -r
