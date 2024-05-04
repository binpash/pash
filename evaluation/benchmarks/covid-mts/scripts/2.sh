#!/bin/bash
# Days a vehicle is on the road

cat "$IN" |                           # assumes saved input
  sed 's/T..:..:..//' |             # hide times
  cut -d ',' -f 3,1 |               # keep only day and bus ID
  sort -u |                         # removing duplicate day-buses
  cut -d ',' -f 2 |                 # keep only bus ID
  sort |                            # preparing for uniq
  uniq -c |                         # count unique dates
  sort -k1n |                       # sort in reverse numerical order
  awk -v OFS="\t" "{print \$2,\$1}" # print first date, then count
