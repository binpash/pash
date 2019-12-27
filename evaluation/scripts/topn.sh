#!/bin/bash

# Identify top 1000 terms in an input

N=1000
INPUT=./input/i1G.txt
cat $INPUT | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${N}q

