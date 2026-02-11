#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

if [ ! -f ./in.csv ]; then
    curl -f 'https://atlas-group.cs.brown.edu/data/covid-mts/in_full.csv.gz'> in.csv.gz
    gzip -d in.csv.gz
fi

if [ ! -f ./in_small.csv ]; then
    curl -f 'https://atlas-group.cs.brown.edu/data/covid-mts/in_small.csv.gz' > in_small.csv.gz
    gzip -d in_small.csv.gz
fi

if [ ! -f ./in_tiny.csv ]; then
    curl -f 'https://atlas-group.cs.brown.edu/data/covid-mts/in_tiny.csv.gz' > in_tiny.csv.gz
    gzip -d in_tiny.csv.gz
fi
