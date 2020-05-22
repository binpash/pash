#!/bin/bash

seq 2005 2005 |
    sed 's;^;http://ndr.md/data/noaa/;' |
    sed 's;$;/;' |
    xargs -n 1 curl -s |
    grep gz |
    tr -s ' ' |
    cut -d ' ' -f9 |
    sed 's;^;http://ndr.md/data/noaa/2005/;'

