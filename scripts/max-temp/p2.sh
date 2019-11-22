#!/bin/bash

cat p1.txt |
  grep gz |
# head -n 1 | # <-- remove this live for full scale
  tr -s ' ' |
  cut -d ' ' -f9 > p2.txt


