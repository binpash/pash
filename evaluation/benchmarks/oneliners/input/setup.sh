#!/bin/bash

set -e

# call the script with its absolute name
cd $(dirname $0)

curl 'http://ndr.md/data/dummy/1G.txt' > 1G.txt
curl 'http://ndr.md/data/dummy/1M.txt' > 1M.txt
# download wamerican-insane dictionary and sort according to machine
curl 'http://ndr.md/data/dummy/dict.txt' | sort > dict.txt

if [ "$#" -eq 1 ] && [ "$1" = "--full" ]; then
  echo Generting full-size inputs
  # FIXME PR: Do we need all of them?

  touch 3G.txt
  for (( i = 0; i < 3; i++ )); do
    cat 1G.txt >> 3G.txt
  done

  touch 10G.txt
  for (( i = 0; i < 10; i++ )); do
    cat 1G.txt >> 10G.txt
  done

  touch 100G.txt
  for (( i = 0; i < 10; i++ )); do
    cat 10G.txt >> 100G.txt
  done
fi
