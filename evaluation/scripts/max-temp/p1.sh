#!/bin/bash

seq 2005 2005 |
    sed 's;^;http://ndr.md/noaa/;' |
    sed 's;$;/;' |
    xargs -n 1 curl -s



