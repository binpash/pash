#!/bin/bash
# Vehicles on the road per day

cat "$IN" |                         # 0: assumes saved input
  sed 's/T..:..:..//' |             # 1/4: hide times
  cut -d ',' -f 1,3 |               # 1/4:keep only day and bus no
  sort -u |                         # 2: remove duplicate records due to time 
  cut -d ',' -f 1 |                 # 3: keep all dates
  sort |                            # 5 sort -m uniq awk: preparing for uniq
  uniq -c |                         # count unique dates
  awk "{print \$2,\$1}" > ${OUT}stdout.txt
