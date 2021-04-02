#!/bin/bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

if [[ "$1" == "-c" ]]; then
    rm -f *.bz2 'in.csv'
    exit
fi

if [ ! -f ./in.csv ]; then
  # yesterday=$(date --date='1 days ago' +'%y-%m-%d')
  # curl https://www.balab.aueb.gr/~dds/oasa-$yesterday.bz2 |
  curl -sf 'https://www.balab.aueb.gr/~dds/oasa-2021-01-08.bz2' | bzip2 -d > in.csv
  if [ $? -ne 0 ]; then
    echo "oasa-2021-01-08.bz2 / bzip2 not available, contact the pash authors"
  fi
fi
