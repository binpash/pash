#!/bin/bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

if [[ "$1" == "-c" ]]; then
    rm -f *.bz2 'in.csv' 'in_small.csv'
    exit
fi

setup_dataset() {
  if [ ! -f ./in.csv ] && [ "$1" = "--full" ]; then
    # yesterday=$(date --date='1 days ago' +'%y-%m-%d')
    # curl https://www.balab.aueb.gr/~dds/oasa-$yesterday.bz2 |
    curl -sf 'https://www.balab.aueb.gr/~dds/oasa-2021-01-08.bz2' | bzip2 -d > in.csv
    if [ $? -ne 0 ]; then
      echo "oasa-2021-01-08.bz2 / bzip2 not available, contact the pash authors"
      exit 1
    fi
  elif [ ! -f ./in_small.csv ] && [ "$1" = "--small" ]; then
    if [ ! -f ./in_small.csv ]; then                                                       
      echo "Generating small-size inputs"                                                  
      # FIXME PR: Do we need all of them?                                                  
      curl -sf 'http://pac-n4.csail.mit.edu:81/pash_data/small/in_small.csv' > in_small.csv
    fi                                                                                     
  fi
}

source_var() {
  if [[ "$1" == "--small" ]]; then
    export IN="input/in_small.csv"
  else
    export IN="input/in.csv"
  fi    
}
