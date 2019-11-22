#!/bin/bash

cat p3.txt |
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > p4.txt
