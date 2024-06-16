#!/bin/bash

# 4.6: piece used the most by Belle
# cat $IN | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

# cat $IN | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3


cat $IN | tail -n 1 >${OUT}stdout.txt