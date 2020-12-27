#!/bin/bash

DIR="scripts/input"
FILE="$DIR/10M.txt"
DICT="$DIR/sorted_dict.txt"

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT -
