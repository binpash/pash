#!/bin/bash
# This script is part of a study on OASA's Telematics
# Diomidis Spinellis and Eleftheria Tsaliki
# https://insidestory.gr/article/noymera-leoforeia-athinas

# # Hours each vehicle is on the road
# <in.csv sed 's/T\(..\):..:../,\1/' |
# awk -F, '!seen[$1 $2 $4] {onroad[$4]++; seen[$1 $2 $4] = 1}
#    END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
# sort -k2n > out1

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
cat in.csv |                    # assumes saved input
  sed 's/T\(..\):..:../,\1/' |  # keep times only
  cut -d ',' -f 1,2,4 |         # keep only time date and bus id
  sort -u |                     # removing duplicate entries
  cut -d ',' -f 3 |             # keep only bus ID
  sort |                        # preparing for uniq
  uniq -c |                     # count hours per bus
  sort -k1n |                   # sort in reverse numerical order
  awk -v OFS="\t" "{print \$2,\$1}"     # print first date, then count

# diff out{1,}
