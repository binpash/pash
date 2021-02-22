#!/bin/bash
# Top-N (1000) terms
# from https://dl.acm.org/doi/10.1145/5948.315654

IN=${IN:-$PASH_TOP/evaluation/benchmarks/expert-oneliners/10G.txt}

cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 100q

