#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

if [ ! -f ./in.csv ]; then
    curl -f 'https://atlas-group.cs.brown.edu/data/covid-mts/in.csv.gz'> in.csv.gz
    gzip -d in.csv.gz
fi

if [ ! -f ./in_small.csv ]; then
    curl -f 'https://atlas-group.cs.brown.edu/data/covid-mts/in_small.csv.gz' > in_small.csv.gz
    gzip -d in_small.csv.gz
fi

if [ ! -f ./in_200M.csv  ]; then
    touch in_200M.csv

    for (( i = 0; i < 25; i++ )); do
        cat in_small.csv >>in_200M.csv
    done

    "$PASH_TOP/scripts/append_nl_if_not.sh" in_200M.csv
fi

if [ ! -f ./in_500M.csv  ]; then
    touch in_500M.csv

    for (( i = 0; i < 62; i++ )); do
        cat in_small.csv >>in_500M.csv
    done

    "$PASH_TOP/scripts/append_nl_if_not.sh" in_500M.csv
fi

if [ ! -f ./in_1G.csv  ]; then
    touch in_1G.csv

    for (( i = 0; i < 125; i++ )); do
        cat in_small.csv >>in_1G.csv
    done

    "$PASH_TOP/scripts/append_nl_if_not.sh" in_1G.csv
fi