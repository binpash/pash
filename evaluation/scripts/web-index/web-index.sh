#!/bin/bash

# Requires: pandoc, node, python

#read from a local url file
# gen() { head -n 1 ./urls.txt; }
# mkfifo s1 s2 s3
# xargs -n 1 curl -s |
# rm s1 s2 s3

cat ./seed.html | 
# N.b.: this has to be in one line, and works only in GNU grep (for -P)
# grep -Po '\b(([\w-]+://?|www[.])[^\s()<>]+(?:\([\w\d]+\)|([^[:punct:]\s]|/)))'
         # '([a-zA-Z][a-zA-Z0-9]*)://([^ /]+)(/?[^ ]*)|([^ @]+)@([^ @]+)'

SEEN="./seen-urls.txt";
SOURCE="./seed.txt";

mkfifo s1 s2
echo './seed.html' | tee s1

# Download page and clone stream as s1, s2
tail -f $SOURCE | xargs -n 1 curl -s | tee s1 s2 >/dev/null 

# s1---URL manipulation: Get all URLs, diff with SEEN, and write back to SOURCE and SEEN
cat s1 | ./grep-url.js | comm -3 $SEEN - | tee -a $SOURCE | sort >> $SEEN # maybe could do >>

# s2---NLP manipulation:  get text
cat s2 | pandoc --from html --to plain | tr -cs A-Za-z '\n' | tr A-Z a-z | grep -vwFf stopwords.txt | ./stem-words.js # stem-to-roots
tee 
  >(tee shifted |  tail +2 | paste shifted - | sort | uniq -c | sort -rn >> 2grams.txt) # 2-grams
  >(tee shifted |  tail +3 | paste shifted - | sort | uniq -c | sort -rn >> 3grams.txt) # 3-grams
  >(tee shifted |  tail +4 | paste shifted - | sort | uniq -c | sort -rn >> 4grams.txt) # 4-grams
 | sort | uniq -c | sort -rn >> 1grams.txt # 1-gram frequencies

# lynx -dump -stdin
rm s1 s2
