#!/bin/bash
# Find all N-Grams (for N==2, bigrams)
alias gen='head -n 2 ./input/input.txt'
# alias gen='cat ./input/input.txt'
mkfifo s2
gen | tr -cs A-Za-z '\n' | tr A-Z a-z | tee s2 |  tail +2 | paste s2 - | sort | uniq > results


