#!/bin/bash

set -e

rm *.txt

echo 1 10 11 12 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 | xargs -n 1 |
  sed 's/$/.txt/' | sed 's;^;http://ndr.md/data/unix50/;' | xargs -n1 wget

if [ "$#" -eq 1 ] && [ "$1" = "--fetch-full" ]; then
  echo 1 10 11 12 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 |
    xargs -n 1 |
    sed 's/$/.1G.txt/' |
    sed 's;^;http://ndr.md/data/unix50/;' |
    xargs -n1 wget
fi

if [ "$#" -eq 1 ] && [ "$1" = "--gen-full" ]; then
  echo Generting full-size inputs

  for file in *.txt; do
    new_file=$(basename $file .txt).1G.txt
    max=$(echo "1000000000 / $(stat --printf="%s\n" $file)" | bc)
    echo "generating 1-G $new_file (${max}x increase)"
    for (( i = 0; i < max ; i++ )); do
      cat $file >> $new_file;
    done
  done
fi
