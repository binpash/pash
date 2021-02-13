#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# Vehicles on the road per day
<input.csv sed 's/T..:..:..//' |
awk -F, '!seen[$1 $3] {onroad[$1]++; seen[$1 $3] = 1}
   END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
sort > 1a.txt

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |              # decompress
# Replace the line below with the two lines above to stream the latest file
cat input.csv |             # assumes saved input
  sed 's/T..:..:..//' |     # hide times
  awk -F, '{print $1,$3}' | # keep only day and bus no
  sort |                    # 
  uniq -c |                 # remove duplicate records due to time
  awk '{print $2}' |        # keep all dates
  sort |
  uniq -c |                 # count unique dates
  awk '{print $2,$1}' |     # print first date, then count
  tr ' ' '\t' > 1b.txt  # replace space w/ tab as per original script

diff 1{a,b}.txt
