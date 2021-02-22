#!/bin/bash
# Sort input

IN=${IN:-$PASH_TOP/evaluation/benchmarks/expert-oneliners/10G.txt}

cat $IN | sort

