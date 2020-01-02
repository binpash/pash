#!/bin/bash

# Requires avoiding buffering lines
# http://mywiki.wooledge.org/BashFAQ/009

# Requires: pandoc, node, python

# iconv is stateless
# ```bash
# echo 'hello there, would you like some umlauts?\n\n ü ä ö ß or some unicode? ☺' |
# iconv -c -t ascii//TRANSLIT
# ```

cat  >seed.txt <<EOM
# Do not start this experiment without redirecting **all** HTTP 
# requests to one of our servers.
https://en.wikipedia.org
EOM

PROXY=$([ "$(hostname)" == "deathstar" ] && echo "gamma.ndr.md" || echo "localhost")
SEEN="./seen.txt";
SOURCE="./seed.txt";
i=0

rm s1 s2
mkfifo s1 s2

echo $PROXY
# Download page and clone stream as s1, s2
# curl http://example.com/2.txt --connect-to '::158.130.4.114:8080'
tail -f $SOURCE |
grep --line-buffered -v '^#' |
sed -u 's/^https/http/' |
head -n 2 |
xargs -n 1 curl -s --connect-to "::${PROXY}:8080"  |
# {
#   # this lambda is only for hitting multiple servers to increase throughput
#   read url;
#   hosts=(localhost) # Can add more servers (gamma.ndr.md delta.ndr.md)
#   sleep 1;
#   echo $url >&2;
#   curl -s  --connect-to "::${hosts[$((i++ % 2))]}:8080" 
# } |
# grep -in dish | # To test if it hits our local HTTP server
# head -n 1 |
  tee s1 >/dev/null &

echo "https://en.wikipedia.org" >> $SOURCE &

# # s1---URL manipulation: Get all URLs, diff with SEEN, and write back to SOURCE and SEEN
# FIXME: sort won't work here because it's gonna wait indefinitely---but { read v; echo $v | comm .. } should work
cat s1 |
./grep-url.js |
# grep -Eoi '<a [^>]+>' | 
# grep -Eo 'href="[^\"]+"' | 
# grep -Eo '(http|https)://[^/"]+' |
grep '^http' |
# tr 'x' 'x' # | tee -a $SOURCE >> $SEEN # maybe could do >>
# {
#   read v;
#   echo "-->  $v" >&2;
#   echo $v; | comm -23 <(cat $SEEN | sort) -;
# } |
grep -v '^#' |
tr 'x' 'x' | tee -a $SOURCE >> $SEEN


# 
# # s2---NLP manipulation:  get text
# cat s2 | pandoc --from html --to plain | tr -cs A-Za-z '\n' | tr A-Z a-z | grep -vwFf stopwords.txt | ./stem-words.js # stem-to-roots
# tee 
#   >(tee shifted |  tail +2 | paste shifted - | sort | uniq -c | sort -rn >> 2grams.txt) # 2-grams
#   >(tee shifted |  tail +3 | paste shifted - | sort | uniq -c | sort -rn >> 3grams.txt) # 3-grams
#   >(tee shifted |  tail +4 | paste shifted - | sort | uniq -c | sort -rn >> 4grams.txt) # 4-grams
#  | sort | uniq -c | sort -rn >> 1grams.txt # 1-gram frequencies
# 
# # lynx -dump -stdin
# rm s1 s2
