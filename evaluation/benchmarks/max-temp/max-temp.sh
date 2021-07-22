#!/bin/bash

FROM=${FROM:-2015}
TO=${TO:-2015}
IN=${IN:-'http://ndr.md/data/noaa/'}
fetch=${fetch:-"curl -s"}

seq $FROM $TO |
  sed "s;^;$IN;" |
  sed 's;$;/;' |
  xargs -r -n 1 $fetch |
  grep gz |
  tr -s ' \n' |
  cut -d ' ' -f9 |
  sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
  sed "s;^;$IN;" |
  xargs -n1 curl -s |
  gunzip |
  ## Processing
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1
