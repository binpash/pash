#!/bin/bash
# Calculate sort twice

IN=${IN:-$PASH_TOP/evaluation/benchmarks/expert-oneliners/10G.txt}

cat $IN | tr A-Z a-z | sort | sort -r
