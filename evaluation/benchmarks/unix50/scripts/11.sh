#!/bin/bash

# 4.5: 4.4 + pawns
cat $IN | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr
