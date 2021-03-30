#!/bin/bash

# another solution for capturing HTTP status code
# https://superuser.com/a/590170

curl -f 'http://ndr.md/data/dummy/1M.txt' > 1M.txt

if [[ $? -ne 0 ]]; then
  curl -f 'http://www.gutenberg.org/files/2600/2600-0.txt' | head -c 1M > 1M.txt
  if [[ $? -ne 0 ]]; then
    echo 'cannot find 1M.txt -- please contact the developers of pash'
    exit 1
  fi
fi

curl -f 'http://ndr.md/data/dummy/all_cmds.txt' > all_cmds.txt
if [[ $? -ne 0 ]]; then
  # This should be OK for tests, no need for abort
  ls /usr/bin/* > all_cmds.txt
fi

curl -f 'http://ndr.md/data/dummy/words' > words

if [[ $(uname) == 'Darwin' ]]; then
  # On OSX, for development
  cp /usr/share/dict/web2 words
  if [[ $? -ne 0 ]]; then
    echo 'cannot find dict file -- please contact the pash developers'
    exit 1
  fi
else
  # On Linux and Debian, for experiments
  # apt install wamerican-insane
  cp /usr/share/dict/words words
  if [[ $? -ne 0 ]]; then
    echo 'cannot find dict file -- please contact the pash developers'
    exit 1
  fi
fi

rm -f 10M.txt
touch 10M.txt
for (( i = 0; i < 10; i++ )); do
  cat 1M.txt >> 10M.txt
done

## Re-sort words for this machine
sort words > sorted_words
