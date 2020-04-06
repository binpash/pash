#!/bin/bash

# A somewhat suboptimal way of calculating 3-grams.
# Part of the intention is to highlight overheads of tagging each stream element

IN=./input/1G.txt
OUT=./output/out.txt

mkfifo s2 s3

cat $IN |
# head -n 2 |
  sed 's/[^a-zA-Z0-9]/ /g' |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  tee s2 |
  tail +2 |
  paste s2 - |            # At this point the stream has two elements <a, b>
  tee s3 |
  cut -f 1 |
  tail +3 |
  paste s3 - |            # Joining (1) the first two words <a, b>, (2) <c>
  sort |
  uniq > $OUT
rm s2 s3



