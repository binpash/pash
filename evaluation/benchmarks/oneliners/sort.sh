#!/bin/bash
# Sort input

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

cat $IN | sort

