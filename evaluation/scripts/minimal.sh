#!/bin/bash
# Timings on my machine (nv):
# 1M: 1.25s
# 1G: > half an hour
IN=./input/i1G.txt # Change G to M for small input
cat $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > ./output/1.txt
