#!/bin/bash

# Word frequencies:

INPUT=./input/i1G.txt
cat ./input.txt | tr -cs A-Za-z'\n' | tr A-Z a-z | sort | uniq -c | sort


