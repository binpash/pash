#!/bin/bash

seq 2005 2005 |
  sed 's;^;http://ndr.md/noaa/;' |
  sed 's;$;/;' |
  xargs -n 1 curl -s > p1.txt

cat p1.txt |
  grep gz |
  head -n 1 | # <-- remove this live for full scale
  tr -s ' ' |
  cut -d ' ' -f9 > p2.txt

cat p2.txt | 
  sed 's;^;http://ndr.md/noaa/2005/;' |
  xargs -n1 curl -s |
  gunzip > p3.txt

cat p3.txt |
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > p4.txt
