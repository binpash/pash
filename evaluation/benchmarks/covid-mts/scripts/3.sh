#!/bin/bash
# Hours each vehicle is on the road

cat "$IN" |                     # assumes saved input
  sed 's/T\(..\):..:../,\1/' |  # keep times only
  cut -d ',' -f 1,2,4 |         # keep only time date and bus id
  sort -u |                     # removing duplicate entries
  cut -d ',' -f 3 |             # keep only bus ID
  sort |                        # preparing for uniq
  uniq -c |                     # count hours per bus
  sort -k 1 -n |                # sort in reverse numerical order
  awk "{print \$2,\$1}"         # print first date, then count
  > "$OUT"stdout.txt                      # save output
