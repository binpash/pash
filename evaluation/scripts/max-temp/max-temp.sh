#!/bin/bash

#NOTE: The `head -n 1 below is for minimizing the number of pages to be seen

# This is similar to {1995..2005}, but this requires shell expansion rules that
# are quite convoluted
seq 2005 2005 | 
  sed 's;^;http://ndr.md/noaa/;' |
  sed 's;$;/;' |
  xargs -n 1 curl -s |
  grep gz |
  head -n 1 | # <-- remove this live for full scale
  tr -s ' ' |
  cut -d ' ' -f9 |
  sed 's;^;http://ndr.md/noaa/2005/;' |
  xargs -n1 curl -s |
  gunzip |
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1

# what can you do with a named function?

# can you use command substitution?

