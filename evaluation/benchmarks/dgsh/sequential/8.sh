#!/bin/bash

# Consistent sorting across machines
export LC_ALL=C

# Convert input into a ranked frequency list
ranked_frequency()
{
	awk '{count[$1]++} END {for (i in count) print count[i], i}' |
	# We want the standard sort here
	sort -rn
}

# Convert standard input to a ranked frequency list of specified n-grams
ngram()
{
	local N=$1

	perl -ne 'for ($i = 0; $i < length($_) - '$N'; $i++) {
		print substr($_, $i, '$N'), "\n";
	}' |
	ranked_frequency
}

# Temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)

cat > "$file1"

# Split input one word per line
tr -cs a-zA-Z '\n' < "$file1" > "$file2"

# Digram frequency
echo "Digram frequency"
ngram 2 < "$file2" 

# Trigram frequency
echo "Trigram frequency"
ngram 3 < "$file2" 

# Word frequency
echo "Word frequency"
ranked_frequency < "$file2"

# Store number of characters to use in awk below
nchars=$(wc -c < "$file1")

# Character frequency
echo "Character frequency"
sed 's/./&\
/g' < "$file1" | 
# Print absolute
ranked_frequency | tee "$file3"

# Print relative
echo "Relative character frequency"
awk -v NCHARS=$nchars 'BEGIN {
		OFMT = "%.2g%%"}
		{print $1, $2, $1 / NCHARS * 100}' "$file3"

