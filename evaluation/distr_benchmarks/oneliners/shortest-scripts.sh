#!/bin/bash
# Find the shortest scripts 
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7
# +p.95 multiple sed
# +p.XX crawler

# FIX: Input here should be a set of commands, more precisely, the ones on this specific machine.

IN=${IN:-/oneliners/all_cmdsx100.txt}

hdfs dfs -cat $IN | xargs file | grep "shell script" | cut -d: -f1 | xargs -L 1 wc -l | grep -v '^0$' | sort -n | head -15
