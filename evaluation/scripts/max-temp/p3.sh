#!/bin/bash

sed 's;^;http://ndr.md/noaa/2005/;' |
    xargs -n1 curl -s |
    gunzip

