#!/bin/bash

## This is now included in the directory by default
# curl 'http://ndr.md/corpus/dummy/1M.txt' > 1M.txt

rm -f 10M.txt
touch 10M.txt
for (( i = 0; i < 10; i++ )); do
  cat 1M.txt >> 10M.txt
done

rm -f 100M.txt
touch 100M.txt
for (( i = 0; i < 10; i++ )); do
  cat 10M.txt >> 100M.txt
done

rm -f 1G.txt
touch 1G.txt
for (( i = 0; i < 10; i++ )); do
  cat 100M.txt >> 1G.txt
done

## Re-sort words for this machine
sort words > sorted_words