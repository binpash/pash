#!/bin/bash

FROM=${FROM:-2015}
TO=${TO:-2015}
IN='https://atlas.cs.brown.edu/data/noaa/'
fetch=${fetch:-"curl -s"}

data_file=temperatures.txt

## Downloading and extracting
seq $FROM $TO |
  sed "s;^;$IN;" |
  sed 's;$;/;' |
  xargs -r -n 1 curl -k -s |
  grep gz |
  head -n 10 |
  tr -s ' \n' |
  cut -d ' ' -f9 |
  sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
  sed "s;^;$IN;" |
  xargs -r -n 1 curl -k -s |
  gunzip > "${data_file}"

## Processing
cat "${data_file}" |
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > max.txt

cat "${data_file}" |
  cut -c 89-92 |
  grep -v 999 |
  sort -n |
  head -n1 > min.txt

cat "${data_file}" |
  cut -c 89-92 |
  grep -v 999 |
  awk "{ total += \$1; count++ } END { print total/count }" > average.txt 
