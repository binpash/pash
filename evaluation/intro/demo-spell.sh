#!/bin/bash

# Would this expansion work:
# cd "$(dirname "${BASH_SOURCE[0]}")"
FILE="$PASH_TOP/evaluation/scripts/input/100M.txt"
DICT="$PASH_TOP/evaluation/scripts/input/sorted_words"

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT -
