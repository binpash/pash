#!/bin/bash

if [ "$#" != "1" ]; then
    echo "Usage: $0 [hwXX]"
    exit 1
fi

if [ -d "$1" ]; then
    echo "Grading directory already exists"
    exit 2
fi

mkdir $1
mkdir $1/submissions
cp ../dropbox/$1/* $1/submissions
