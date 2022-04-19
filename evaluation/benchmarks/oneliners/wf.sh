#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

IN=$PASH_TOP/evaluation/benchmarks/oneliners/input/1M.txt

cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn 
