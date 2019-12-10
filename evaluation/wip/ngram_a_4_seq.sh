#!/bin/bash
# Find all N-Grams (for N==2, bigrams)
cat $IN $IN $IN $IN | tr -cs A-Za-z '\n' | tr A-Z a-z
