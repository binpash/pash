#!/bin/bash

# A bash script for determining how many pages are in a folder of OpenOffice documents
# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7

# Require: libimage-exiftool-perl, bc
# Data:
#   http://ndr.md/data/dummy/large.pdf
# More data:
#   https://arxiv.org/help/bulk_data

IN=./input/large.pdf
OUT=./output/out.txt

echo "$(exiftool $IN |
    grep Page-count |
    cut -d ":" -f2 |
    tr '\n' '+')""0" |
  bc |
  sed 's/^/\n/' > $OUT
