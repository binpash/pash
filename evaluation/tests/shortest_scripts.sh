#!/bin/bash
# A bash script for finding the shortest scripts 
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7
# +p.95 multiple sed
# +p.XX crawler

# cut -d: -f1 -> cut -d : -f 1; as parser recognizes option arguments only if given with whitespace
# head -15 -> head -n 15; not documented in man page 
cat $IN | xargs file | grep "shell script" | cut -d : -f 1 | xargs -L 1 wc -l | grep -v '^0$' | sort -n | head -n 15
