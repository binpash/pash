#!/bin/bash

# Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)

# Save the ls output to a temporary file
ls -n > "$file1"

# Reorder fields in DIR-like way
awk '!/^total/ {print $6, $7, $8, $1, sprintf("%8d", $5), $9}' "$file1" > "$file2"

# Count number of files
wc -l "$file1" | tr -d \\n > "$file3"
echo -n ' File(s) ' >> "$file3"
awk '{s += $5} END {printf("%d bytes\n", s)}' "$file1" >> "$file3"

# Count number of directories and print label for number of dirs and calculate free bytes
grep -c '^d' "$file1" | tr -d \\n > "$file4"
df -h . | awk '!/Use%/{print " Dir(s) " $4 " bytes free"}' >> "$file4"

# Display the results
cat "$file2" "$file3" "$file4"
