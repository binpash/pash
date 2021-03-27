#!/bin/bash

# https://gist.github.com/noamross/86fba413e0769069e3955d1c9bc530ae
funzip $1|  # uncompress first file in zip
tr -d '\000' | #remove null characters
sed "/^\s*$/d; s/ \{1,\}\t/\t/g; s/\t \{1,\}/\t/g; s/\r//" |  #removes empty lines, whitespace around tabs, extra newlines
cut -s -f 1,3,4,5,6,8,12,13,14,15,16,17,18,19,20,21,23,24,25,26,34,35,36,38,40,42,44,45,46,85,86,87,88,89 #| #only select certain columns
pv -N Process -c  |
gzip -9 | 
pv -N Compress -c > $1.gz
