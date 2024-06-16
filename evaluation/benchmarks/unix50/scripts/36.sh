#!/bin/bash

# 11.2: most repeated first name in the list?
# cat $IN | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w 1 | sed 1d

cat $IN | fmt -w 1 >${OUT}stdout.txt

