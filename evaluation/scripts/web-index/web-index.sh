#!/bin/bash

# Requires avoiding buffering lines
# http://mywiki.wooledge.org/BashFAQ/009

# Requires: pandoc, node, python
# run `npm install` to install packages

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

rm s1 s2 drainbuf
mkfifo s1 s2 drainbuf
echo $PROXY

# TODO (nv): I have a feeling this creates duplicates
cat drainbuf |
  tee /dev/null >&2 & 

# to enable line buffering: stdbuf -oL
tail -f $SOURCE |
# removes buffering
tee drainbuf |
grep --line-buffered -v '^#' |
sed -u 's/^https/http/' |
head -n 2 |
# to debug curl, remove `-s`
xargs -0 -n 1 -d '\n' curl -s --connect-to "::${PROXY}:8080"  |
# {
#   # this lambda is only for hitting multiple servers to increase throughput
#   read url;
#   hosts=(localhost) # Can add more servers (gamma.ndr.md delta.ndr.md)
#   sleep 1;
#   echo $url >&2;
#   curl -s  --connect-to "::${hosts[$((i++ % 2))]}:8080" 
# } |
  tee -a s1 -a s2 >/dev/null &

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
tee -a $SEEN -a $SOURCE > /dev/null &

rm shift{1,2,3} {1,2,3}grams
mkfifo shift{1,2,3} {1,2,3}grams

cat '1grams' |
    sort |
    uniq -c |
    sort -rn >> 1-grams.txt &

cat '2grams' |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    tee shift1 |
    tail +2 |
    paste shift1 - |
    sort |
    uniq -c |
    sort -rn >> 2-grams.txt &

cat '3grams' |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    tee shift2 |
    tail +2 |
    paste shift2 - |
    tee shift3 |
    cut -f 1 |
    tail +3 |
    paste shift3 - |
    sort |
    uniq -c |
    sort -rn >> 3-grams.txt &

# s2---NLP manipulation:  get text
cat s2 |
  pandoc --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  iconv -c -t ascii//TRANSLIT |
  grep -vwFf stopwords.txt |
  ./stem-words.js | # stem-to-roots
  tee '3grams' '2grams' '1grams' > /dev/null &

wait
# lynx -dump -stdin
# rm s1 s2
