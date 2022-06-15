#!/bin/bash
FROM=${FROM:-2015}
TO=${TO:-2015}
IN=${IN:-'http://ndr.md/data/noaa/'}
fetch=${fetch:-"curl -s"}

data_file=temperatures.txt

if [[ "$1" == "--extended" ]]; then
  echo "Downloading extended input"
  dataset_size=14418
else
  dataset_size=1442
fi

## Downloading and extracting
seq $FROM $TO |
  sed "s;^;$IN;" |
  sed 's;$;/;' |
  xargs -r -n 1 $fetch |
  grep gz |
  tr -s ' \n' |
  cut -d ' ' -f9 |
  sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
  sed "s;^;$IN;" |
  head -n $dataset_size |
  xargs -n1 $fetch |
  gunzip > "${data_file}"

hdfs dfs -mkdir /max-temp
hdfs dfs -put "${data_file}" /max-temp/"${data_file}"