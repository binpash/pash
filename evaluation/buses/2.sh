#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# Days a vehicle is on the road
<input.csv sed 's/T..:..:..//' |
awk -F, '!seen[$1 $3] {onroad[$3]++; seen[$1 $3] = 1}
   END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
sort -k2n >2a.txt

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
cat input.csv |                 # assumes saved input
  sed 's/T..:..:..//' |         # hide times
  awk -F, '{print $3,$1}' |     # keep only day and bus no
  sort -n |                     # 
  uniq -c |                     # remove duplicate records due to time
  awk '{print $2}' |            # keep all buses
  sort |
  uniq -c |                     # count unique dates
  awk '{print $2,$1}' |         # print first date, then count
  sort -k2n |                   # sort in reverse numerical order
  tr ' ' '\t' > 2b.txt      # replace space w/ tab as per original script

diff 2{a,b}.txt
