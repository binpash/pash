#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
. "$PASH_TOP/scripts/utils.sh"
cd $(dirname $0)

[ "$1" = "-c" ] && rm-files 100M.txt words sorted_words


if [ ! -f ./100M.txt ]; then
  curl -sf --connect-timeout 10 'ndr.md/data/dummy/100M.txt' > 100M.txt
  if [ $? -ne 0 ]; then
    # Pipe curl through tac (twice) in order to consume all the output from curl.
    # This way, curl can write the whole page and not emit an error code.
    curl -fL 'http://www.gutenberg.org/files/2600/2600-0.txt' | tac | tac | head -c 1M > 1M.txt
    [ $? -ne 0 ] && eexit 'cannot find 1M.txt'
    touch 100M.txt
    for (( i = 0; i < 100; i++ )); do
      cat 1M.txt >> 100M.txt
    done
  fi
  append_nl_if_not ./100M.txt
fi

if [ ! -f ./words ]; then
  curl -sf --connect-timeout 10 'http://ndr.md/data/dummy/words' > words
  if [ $? -ne 0 ]; then
    curl -sf 'https://zenodo.org/record/7650885/files/words' > words
    if [ $? -ne 0 ]; then
      if [ $(uname) = 'Darwin' ]; then
        cp /usr/share/dict/web2 words || eexit "cannot find dict file"
      else
        # apt install wamerican-insane
        cp /usr/share/dict/words words || eexit "cannot find dict file"
      fi
    fi
  fi
  append_nl_if_not words
fi

## Re-sort words for this machine
if [ ! -f ./sorted_words ]; then
  sort words > sorted_words
fi
