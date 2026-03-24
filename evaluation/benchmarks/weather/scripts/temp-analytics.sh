#!/bin/bash

## Processing
cat $IN |
  cut -c 88-92 |
  grep -v 999 |
  sort -rn |
  head -n 1 > ${OUT}max.stdout.txt

cat $IN |
  cut -c 88-92 |
  grep -v 999 |
  sort -n |
  head -n 1 > ${OUT}min.stdout.txt

cat $IN |
  cut -c 88-92 |
  grep -v 999 |
  awk "{ total += \$1; count++ } END { print total/count }" > ${OUT}average.stdout.txt 