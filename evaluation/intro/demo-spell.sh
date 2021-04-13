#!/bin/bash

# Would this expansion work:
# cd "$(dirname "${BASH_SOURCE[0]}")"
FILE="$PWD/input/100M.txt"
DICT="$PWD/input/sorted_words"

echo $FILE $DICT

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT -
