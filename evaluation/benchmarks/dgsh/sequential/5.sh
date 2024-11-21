#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)

export LC_ALL=C

cat > "$file1"

# Find errors

# Obtain list of words in text
cat "$file1" |
tr -cs A-Za-z \\n |
tr A-Z a-z |
sort -u > "$file2"

# Ensure dictionary is compatibly sorted
cat "$file1" |
sort /usr/share/dict/words > "$file3"

# List errors as a set difference
comm -23 "$file2" "$file3" > "$file4"

fgrep -f "$file4" -i --color -w -C 2 "$file1"