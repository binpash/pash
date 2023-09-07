#!/bin/bash

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)
file5=$(mktemp)

cat > $file1

# Consistent sorting across machines
export LC_ALL=C

# Stream input from file and split input one word per line
tr -cs a-zA-Z '\n' < "$file1" |
# Create list of unique words
sort -u > "$file2"

# List two-letter palindromes
sed 's/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t
	g' "$file2" > "$file3"

# List four consecutive consonants
sed -E 's/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t
	g' "$file2" > "$file4"

# List length of words longer than 12 characters
awk '{if (length($1) > 12) print "l:", length($1);
	else print ""}' "$file2" > "$file5"

# Paste the four streams side-by-side
paste "$file2" "$file3" "$file4" "$file5" | 
# List only words satisfying one or more properties
fgrep :
