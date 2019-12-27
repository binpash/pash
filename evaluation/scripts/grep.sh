#!/bin/bash

# Substring search  of  a complex  string, using  backtracking  (a superset  of
# Thompson's high-perf NFA used in GNU)

# Data: ./input/i1G.txt

cat ./input/i1G.txt |
  grep '[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}' 


