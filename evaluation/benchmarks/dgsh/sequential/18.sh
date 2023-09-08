#!/bin/bash

# Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

# Read the input stream and save to a temporary file
cat $INPUT_FILE > "$file1"

# Process the input in two different ways
cut -d , -f 5-6 "$file1" > "$file2"
cut -d , -f 2-4 "$file1" > "$file3"

# Merge the processed results
paste -d , "$file2" "$file3"
