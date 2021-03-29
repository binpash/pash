#!/bin/bash
# Top-N (1000) terms
# from https://dl.acm.org/doi/10.1145/5948.315654

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 100q

