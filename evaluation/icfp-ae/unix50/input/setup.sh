#!/bin/bash

set -e

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

## FIXME: These inputs are already 1G when downloaded
## FIXME: Also, wget is not silent like curl in the other setup scripts.

inputs=(
  1 10 11 12 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9
)
if [[ "$1" == "-c" ]]; then
  for input in ${inputs[@]}
  do
    rm -f "${input}.txt"
  done
  exit
fi

for input in ${inputs[@]}
do
  if [ ! -f "${input}.txt" ]; then
    wget "http://ndr.md/data/unix50/${input}.txt"
    "$PASH_TOP/scripts/append_nl_if_not.sh" "${input}.txt"
  fi
done

## FIXME: Calling this script with --full is not idempotent.
if [ "$#" -eq 1 ] && [ "$1" = "--full" ]; then
  for file in *.txt; do
    echo '' > temp.txt
    for (( i = 0; i < 10; i++ )); do
      cat $file >> temp.txt
    done
    mv temp.txt $file
  done
fi
