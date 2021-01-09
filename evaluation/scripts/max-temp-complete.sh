#!/bin/bash

seq 2000 2004 |
    sed 's;^;http://ndr.md/data/noaa/;' |
    sed 's;$;/;' |
    xargs -n 1 curl -s |
    grep gz |
    tr -s ' ' |
    cut -d ' ' -f9 |
    head |
    sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
    ## TODO: Should this be 2005?
    sed 's;^;http://ndr.md/data/noaa/2005/;' |
    xargs -n1 curl -s |
    gunzip |
    ## Processing
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1