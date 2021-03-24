#!/bin/bash
# Match complex regular-expression over input

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

cat $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
