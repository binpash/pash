#!/bin/bash

if [ "$#" != "1" ]; then
    echo "Usage: $0 [hwXX]"
    exit 1
fi

if [ ! -d "$1/mail" ]; then
    echo "Couldn't find mail directory (looked in $1/grading)"
    exit 2
fi

cd $1/mail

for s in `ls`; do
    ../../mail.scpt "[cs131] $1 grade" $s
done
