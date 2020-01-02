#!/bin/bash

# Requires avoiding buffering lines
# http://mywiki.wooledge.org/BashFAQ/009

# Requires: pandoc, node, python

# iconv is stateless
# ```bash
# echo 'hello there, would you like some umlauts?\n\n ü ä ö ß or some unicode? ☺' |
# iconv -c -t ascii//TRANSLIT
# ```

# https://webcache.googleusercontent.com/search?q=cache:l-R5pHOys9MJ:https://stackoverflow.com/questions/207047/what-linux-unix-commands-are-outdated-and-have-powerful-alternatives+&cd=13&hl=en&ct=clnk&gl=us
# https://github.com/flonatel/pipexec

cat  >seed.txt <<EOM
# Do not start this experiment without redirecting **all** HTTP 
# requests to one of our servers.
http://en.wikipedia.org
EOM

PROXY=$([ "$(hostname)" == "deathstar" ] && echo "gamma.ndr.md" || echo "localhost")
SEEN="./seen.txt";
SOURCE="./seed.txt";
i=0

rm s1 s2
mkfifo s1 s2

echo $PROXY

# to enable line buffering: stdbuf -oL
tail -f $SOURCE |
# removes buffering
tee >(tee /dev/null >&2) |
grep --line-buffered -v '^#' |
sed -u 's/^https/http/' |
head -n 2 |
# to debug curl, remove `-s`
xargs -0 -n 1 -d '\n' curl --connect-to "::${PROXY}:8080"  |
# {
#   # this lambda is only for hitting multiple servers to increase throughput
#   read url;
#   hosts=(localhost) # Can add more servers (gamma.ndr.md delta.ndr.md)
#   sleep 1;
#   echo $url >&2;
#   curl -s  --connect-to "::${hosts[$((i++ % 2))]}:8080" 
# } |
  tee s1 >/dev/null &

echo "http://en.wikipedia.org" >> $SOURCE &

# # s1---URL manipulation: Get all URLs, diff with SEEN, and write back to SOURCE and SEEN
# FIXME: sort won't work here because it's gonna wait indefinitely---but { read v; echo $v | comm .. } should work
cat s1 |
./grep-url.js |
# grep -Eoi '<a [^>]+>' | 
# grep -Eo 'href="[^\"]+"' | 
# grep -Eo '(http|https)://[^/"]+' |
grep '^http' |
# filtering in-place, because we can't do sort
# realized that xargs is a window operator!
# it's impossible to use <( ) with xargs because it's evaluated before 
# xargs expands `{}` --- so the solution is to call bash
# quotes are used to avoid command injection
# TODO: Alternatively, this can be written as a (possibly anon.) function
xargs -0 -n 1 -d '\n' -I {} bash -c 'echo "$(grep -Fxv -f ./test.txt <(echo "{}"))"' |
# grep -w amp
# {
#   read v;
# # echo "-->  $v" >&2;
#   echo $v # |  -
# } |
grep --line-buffered -v '^#' |
sed -u 's/^https/http/' |
# # tr 'x' 'x' |
# Note the append on SEEN and write on SOURCE
tee -a $SEEN -a $SOURCE > /dev/null

#  # s2---NLP manipulation:  get text
#  cat s2 | pandoc --from html --to plain | tr -cs A-Za-z '\n' | tr A-Z a-z | grep -vwFf stopwords.txt | ./stem-words.js # stem-to-roots
#  tee 
#    >(tee shifted |  tail +2 | paste shifted - | sort | uniq -c | sort -rn >> 2grams.txt) # 2-grams
#    >(tee shifted |  tail +3 | paste shifted - | sort | uniq -c | sort -rn >> 3grams.txt) # 3-grams
#    >(tee shifted |  tail +4 | paste shifted - | sort | uniq -c | sort -rn >> 4grams.txt) # 4-grams
#   | sort | uniq -c | sort -rn >> 1grams.txt # 1-gram frequencies
#  
#  # lynx -dump -stdin
#  rm s1 s2
