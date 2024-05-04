#!/bin/bash
# Vehicles on the road per day

cat "$IN" |                           # assumes saved input
  sed 's/T..:..:..//' |             # hide times
  cut -d ',' -f 1,3 |               # keep only day and bus no
  sort -u |                         # remove duplicate records due to time
  cut -d ',' -f 1 |                 # keep all dates
  sort |                            # preparing for uniq
  uniq -c |                         # count unique dates
  awk -v OFS="\t" "{print \$2,\$1}" # print first date, then count
