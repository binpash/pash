#!/bin/bash

cd "$(dirname "$0")" 

FROM=${FROM:-2015}
TO=${TO:-2015}
IN=${IN:-'https://atlas-group.cs.brown.edu/data/noaa/'}
fetch=${fetch:-"curl -s"}

data_dir=./inputs

## Downloading and extracting
for year in $(seq $FROM $TO)
do
    url_year="$IN$year/"
    data_file="${data_dir}/temperatures.$year.txt"
    ## Note: I am concerned about the use of -s because it might hide some error which might make it hard to find
    curl -s "$url_year" |
    grep gz |
    ## --------- NECESSARY WHEN DOWNLOADING FROM NOAA -----------
    ## note: regexp generated with ChatGPT
    # sed -n 's/.*href="\([^"]*\.gz\)".*/\1/p' |
    # sed "s;^;$url_year;" |
    ## ----------------------------------------------------------
    tr -s ' \n' |
    cut -d ' ' -f9 |
    sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
    sed  "s;^;$IN;" |
    ## Note: I am concerned about the use of -s because it might hide some error which might make it hard to find
    xargs -n1 curl -s |
    gunzip > "${data_file}"
done