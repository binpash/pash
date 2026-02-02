#!/bin/bash

sed 's;^;atlas.cs.brown.edu/data/noaa/;' |
    sed 's;$;/;' |
    xargs -r -n 1 curl -s |
    grep gz |
    tr -s ' \n' |
    cut -d ' ' -f9 |
    sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
    sed 's;^;atlas.cs.brown.edu/data/noaa/;' |
    xargs -n1 curl -s |
    gunzip
