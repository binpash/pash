#!/bin/bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

if [[ "$1" == "-c" ]]; then
    rm -f *.bz2 'in.csv'
    exit
fi

setup_dataset() {
  if [ ! -f ./in.csv ]; then
    # yesterday=$(date --date='1 days ago' +'%y-%m-%d')
    # curl https://www.balab.aueb.gr/~dds/oasa-$yesterday.bz2 |
    curl -sf 'https://www.balab.aueb.gr/~dds/oasa-2021-01-08.bz2' | bzip2 -d > in.csv
    if [ $? -ne 0 ]; then
      echo "oasa-2021-01-08.bz2 / bzip2 not available, contact the pash authors"
      exit 1
    fi
    "$PASH_TOP/scripts/append_nl_if_not.sh"  in.csv
    len=$(cat in.csv | wc -l)
    half_size=$(( $len / 4 ))
    head -n $half_size in.csv > in_small.csv
  fi
}

source_var() {
  if [[ "$1" == "--small" ]]; then
    export IN="input/in_small.csv"
  else
    export IN="input/in.csv"
  fi    
}
