#!/bin/bash

set -e

# call the script with its absolute name
cd $(dirname $0)

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

# another solution for capturing HTTP status code
# https://superuser.com/a/590170

input_files="1M.txt 1G.txt dict.txt 3G.txt 10G.txt 100G.txt"

if [[ "$1" == "-c" ]]; then
  rm -f $input_files
  exit
fi


if [ ! -f ./1M.txt ]; then
  curl -sf 'http://ndr.md/data/dummy/1M.txt' > 1M.txt
  if [ $? -ne 0 ]; then
    echo 'cannot find 1M.txt -- please contact the developers of pash'
    exit 1
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" ./1M.txt
fi

if [ ! -f ./1G.txt ]; then
  curl -sf 'http://ndr.md/data/dummy/1G.txt' > 1G.txt
  if [ $? -ne 0 ]; then
    echo 'cannot find 1G.txt -- please contact the developers of pash'
    exit 1
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" ./1G.txt
fi

# download wamerican-insane dictionary and sort according to machine
if [ ! -f ./dict.txt ]; then
  curl -sf 'http://ndr.md/data/dummy/dict.txt' | sort > dict.txt
  if [ $? -ne 0 ]; then
    echo 'cannot find dict.txt -- please contact the developers of pash'
    exit 1
  fi
  "$PASH_TOP/scripts/append_nl_if_not.sh" ./dict.txt
fi

if [ "$#" -eq 1 ] && [ "$1" = "--full" ]; then
  echo Generting full-size inputs
  # FIXME PR: Do we need all of them?

  if [ ! -f ./3G.txt ]; then
    touch 3G.txt
    for (( i = 0; i < 3; i++ )); do
      cat 1G.txt >> 3G.txt
    done
  fi

  if [ ! -f ./10G.txt ]; then
    touch 10G.txt
    for (( i = 0; i < 10; i++ )); do
      cat 1G.txt >> 10G.txt
    done
  fi

  if [ ! -f ./100G.txt ]; then
    touch 100G.txt
    for (( i = 0; i < 10; i++ )); do
      cat 10G.txt >> 100G.txt
    done
  fi
fi
