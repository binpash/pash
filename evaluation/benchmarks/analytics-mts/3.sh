#!/bin/bash
# Hours each vehicle is on the road

# <in.csv sed 's/T\(..\):..:../,\1/' |
# awk -F, '!seen[$1 $2 $4] {onroad[$4]++; seen[$1 $2 $4] = 1}
#    END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
# sort -k2n > out1

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 | 
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
cat in.csv |                    # assumes saved input
  sed 's/T\(..\):..:../,\1/' |  # keep times only
  awk -F, '{print $1,$2,$4}' |  # keep only day and bus no
  sort |                        # preparing for uniq
  uniq -c |                     # remove duplicate records due to time
  awk '{print $4}' |            # keep all buses
  sort |                        # preparing for uniq
  uniq -c |                     # count unique dates
  awk '{print $2,$1}' |         # print first date, then count
  sort -k2n |                   # sort in reverse numerical order
  tr ' ' '\t' > out             # replace space w/ tab as per original script

# diff out{1,}
