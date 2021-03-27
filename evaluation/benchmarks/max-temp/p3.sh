#!/bin/bash

sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
sed 's;^;http://ndr.md/data/noaa/;'
