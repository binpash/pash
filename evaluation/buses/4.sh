#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# Hours monitored each day
<input.csv sed 's/T\(..\):..:../,\1/' |
awk -F, '!seen[$1 $2] {hours[$1]++; seen[$1 $2] = 1}
   END { OFS = "\t"; for (d in hours) print d, hours[d]}' | 
  sort >4a.txt

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
cat input.csv |                 # assumes saved input
  sed 's/T\(..\):..:../,\1/' |  # keep times only
  awk -F, '{print $1,$2}' |  # keep only day and bus no
  sort |                        # 
  uniq -c |                     # remove duplicate records due to time
  awk '{print $2}' |            # keep all buses
  sort |
  uniq -c |                     # count unique dates
  awk '{print $2,$1}' |         # print first date, then count
  sort |                   # sort in reverse numerical order
  tr ' ' '\t' > 4b.txt      # replace space w/ tab as per original script

diff 4{a,b}.txt
