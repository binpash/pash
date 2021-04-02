#!/bin/bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}


# another solution for capturing HTTP status code
# https://superuser.com/a/590170

if [[ "$1" == "-c" ]]; then
  rm -f 1M.txt all_cmds.txt words sorted_words 10M.txt
  exit
fi

if [ ! -f ./1M.txt ]; then
  curl -sf 'http://ndr.md/data/dummy/1M.txt' > 1M.txt
  if [ $? -ne 0 ]; then
    curl -sf 'http://www.gutenberg.org/files/2600/2600-0.txt' | head -c 1M > 1M.txt
    if [ $? -ne 0 ]; then
      echo 'cannot find 1M.txt -- please contact the developers of pash'
      exit 1
    fi
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" ./1M.txt
fi

if [ ! -f ./all_cmds.txt ]; then
  curl -sf 'http://ndr.md/data/dummy/all_cmds.txt' > all_cmds.txt
  if [ $? -ne 0 ]; then
    # This should be OK for tests, no need for abort
    ls /usr/bin/* > all_cmds.txt
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" ./all_cmds.txt
fi

if [ ! -f ./words ]; then
  curl -sf 'http://ndr.md/data/dummy/words' > words
  if [ $? -ne 0 ]; then
    if [ $(uname) = 'Darwin' ]; then
      # On OSX, for development
      cp /usr/share/dict/web2 words
      if [ $? -ne 0 ]; then
        echo 'cannot find dict file -- please contact the pash developers'
        exit 1
      fi
    else
      # On Linux and Debian, for experiments
      # apt install wamerican-insane
      cp /usr/share/dict/words words
      if [ $? -ne 0 ]; then
        echo 'cannot find dict file -- please contact the pash developers'
        exit 1
      fi
    fi
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" words
fi

if [ ! -f ./10M.txt ]; then
  touch 10M.txt
  for (( i = 0; i < 10; i++ )); do
    cat 1M.txt >> 10M.txt
  done
fi

## Re-sort words for this machine
if [ ! -f ./sorted_words ]; then
  sort words > sorted_words
fi
