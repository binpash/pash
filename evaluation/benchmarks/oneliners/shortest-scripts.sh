#!/bin/bash
# Find the shortest scripts 
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7
# +p.95 multiple sed
# +p.XX crawler

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt}

cat $IN | xargs file | grep "shell script" | cut -d: -f1 | xargs -L 1 wc -l | grep -v '^0$' | sort -n | head -15
