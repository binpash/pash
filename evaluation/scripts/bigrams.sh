#!/bin/bash

# Find all 2-grams in a piece of text

alias gen='head -n 2 ./input/input.txt'
# alias gen='cat ./input/input.txt'
mkfifo s2
gen | tr -cs A-Za-z '\n' | tr A-Z a-z | tee s2 |  tail +2 | paste s2 - | sort | uniq > results


