#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

# Find object files and print defined symbols
find "$1" -name "*.o" | xargs nm > "$file1"

# List all defined (exported) symbols
awk 'NF == 3 && $2 ~ /[A-Z]/ {print $3}' "$file1" | sort > "$file2"

# List all undefined (imported) symbols
awk '$1 == "U" {print $2}' "$file1" | sort > "$file3"

# Print exports that are not imported
comm -23 "$file2" "$file3"
