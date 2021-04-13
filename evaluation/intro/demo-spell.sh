#!/bin/bash

# Would this expansion work:
# cd "$(dirname "${BASH_SOURCE[0]}")"
FILE="$(dirname $0)/input/100M.txt"
DICT="$(dirname $0)/input/sorted_words"

echo $FILE $DICT

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT -
