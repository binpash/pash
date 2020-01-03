#!/bin/bash
# A bash script for finding the shortest scripts 
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7
# +p.95 multiple sed
# +p.XX crawler

cat $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN $IN | xargs file | grep "shell script" | cut -d: -f1 | xargs wc -l | sort -n | head -15
