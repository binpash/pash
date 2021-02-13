#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# Hours each vehicle is on the road
<input.csv sed 's/T\(..\):..:../,\1/' |
awk -F, '!seen[$1 $2 $4] {onroad[$4]++; seen[$1 $2 $4] = 1}
   END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
sort -k2n >3a.txt

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
cat input.csv |           # assumes saved input
  sed 's/T\(..\):..:../,\1/' |  # keep times only
  awk -F, '{print $1,$2,$4}' |  # keep only day and bus no
  sort |                        # 
  uniq -c |                     # remove duplicate records due to time
  awk '{print $4}' |            # keep all buses
  sort |
  uniq -c |                     # count unique dates
  awk '{print $2,$1}' |         # print first date, then count
  sort -k2n |                   # sort in reverse numerical order
  tr ' ' '\t' > 3b.txt      # replace space w/ tab as per original script

diff 3{a,b}.txt
