#!/usr/bin/env dgsh
#
# SYNOPSIS Text properties
# DESCRIPTION
# Read text from the standard input and create files
# containing word, character, digram, and trigram frequencies.
#
# Demonstrates the use of scatter blocks without output and the use
# of stores within the scatter block.
#
# Example:
# curl ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt | text-properties
#
#  Copyright 2013 Diomidis Spinellis
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

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

export -f ranked_frequency
export -f ngram

tee |
{{
	# Split input one word per line
	tr -cs a-zA-Z \\n |
	tee |
	{{
		# Digram frequency
		call 'ngram 2 >digram.txt'
		# Trigram frequency
		call 'ngram 3 >trigram.txt'
		# Word frequency
		call 'ranked_frequency >words.txt'
	}}

	# Store number of characters to use in awk below
	wc -c |
	dgsh-writeval -s nchars

	# Character frequency
	sed 's/./&\
/g' |
	# Print absolute
	call 'ranked_frequency' |
	awk 'BEGIN {
		"dgsh-readval -l -x -q -s nchars" | getline NCHARS
		OFMT = "%.2g%%"}
		{print $1, $2, $1 / NCHARS * 100}' > character.txt
}}
