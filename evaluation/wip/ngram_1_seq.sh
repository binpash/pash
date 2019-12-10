#!/bin/bash
# Find all N-Grams (for N==2, bigrams)
mkfifo s2
cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | tee s2 |  tail +2 | paste s2 - | sort | uniq
