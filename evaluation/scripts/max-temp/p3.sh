#!/bin/bash

cat p2.txt | 
  sed 's;^;http://ndr.md/noaa/2005/;' |
  xargs -n1 curl -s |
  gunzip > p3.txt

