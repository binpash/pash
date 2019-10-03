#!/bin/bash
# Find all N-Grams (for N==2, bigrams)
cat ./input.txt | tr -cs A-Za-z '\n' | tr A-Z a-z > tokens.txt && tail +2 tokens.txt > next.txt && paste tokens.txt next.txt > bigrams.txt && cat bigrams.txt | sort | uniq > results


