#!/bin/bash
# tag: find dublicate files
# from: https://www.spinellis.gr/sw/dgsh/#duplicate-files
#FIXME need dataset 
set -e
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}
rm -f a
mkfifo a
find . -type f |               # Create list of files
  xargs openssl md5 |          # MD5(filename)= 811bfd4b5974f39e986ddc037e1899e7
  sed 's/^MD5(//;s/)= / /' |   # Convert each line into a "filename md5sum" pair
  sort -k2 |                   # Sort by MD5 sum
  tee a |                      # Create a different stream
  awk '{print $2}' |           # Keep only hash
  uniq -d |                    # Keep only duplicates
  join -2 2 - a |              # Join repeated MD5 sums with corresponding file names
  awk ' BEGIN {ORS=""} $0 != prev && prev {print "\n"} END {if (prev) print "\n"} {if (prev) print " "; prev = $1; print $2}'
# Output same files on a single line
rm a
