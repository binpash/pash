#!/bin/bash

seq 2000 2004 |
    sed 's;^;http://ndr.md/data/noaa/;' |
    sed 's;$;/;' |
    xargs -n 1 curl -s



