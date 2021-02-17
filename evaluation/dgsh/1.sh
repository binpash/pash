#!/bin/bash
# Adapted from the DGSH 
# https://www.spinellis.gr/sw/dgsh/#compress-compare

# input wget https://dblp.uni-trier.de/xml/dblp.xml.gz

mkfifo a b c d
cat a | wc -c | sed 's/^/Original size:\t/' &
cat b | xz -c | wc -c | sed 's/^/xz:\t\t/' &
cat c | bzip2 -c | wc -c | sed 's/^/bzip2:\t\t/' &
cat d | gzip -c | wc -c | sed 's/^/gzip:\t\t/' &
cat dblp.xml | tee a b c d >/dev/null
rm a b c d
