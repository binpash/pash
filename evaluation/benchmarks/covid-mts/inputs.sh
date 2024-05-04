#!/bin/bash

cd "$(dirname "$0")" || exit 1

inputs_dir="inputs"

# Function to display usage
usage() {
  echo "Usage: $0 [--small] [--help]"
  echo "  --small      : Download the small dataset"
  echo "  --help       : Show this message and exit"
}

# Function to set up the dataset
setup_dataset() {
  if [[ ! -f "$inputs_dir/in_small.csv" && "$1" == "--small" ]]
  then
    echo "Downloading and decompressing small-size dataset..."
    # TODO: Not implemented; link is broken
    # curl -sf 'http://pac-n4.csail.mit.edu:81/pash_data/small/in_small.csv' > "$inputs_dir/in_small.csv"
  elif [[ ! -f "$inputs_dir/in.csv" ]]
  then
    echo "Downloading and decompressing full-size dataset..."
    curl -sf 'https://www.balab.aueb.gr/~dds/oasa-2021-01-08.bz2' | bzip2 -d > "$inputs_dir/in.csv"
  fi
}

case "$1" in
  --help)
    usage
    ;;
  --small)
    setup_dataset --small
    ;;
  *)
    setup_dataset --full
    ;;
esac

exit 0
