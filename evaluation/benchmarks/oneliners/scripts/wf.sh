#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

cd "$(dirname "$0")" || exit 1

SIZE=500M

IN="input/$SIZE.txt"

cat $IN | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn
