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

IN=${IN:-/usr/share/dict/words}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUTPUT}
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

rm -f a b c d e f g digram.txt trigram.txt words.txt teea_res teeb_res teea teeb
mkfifo a b c d e f g digram.txt trigram.txt words.txt teea_res teeb_res teea teeb

# redirect the input to these 2 processes
# Digram frequency
cat b | ngram 2 |cat > digram.txt &
# Trigram frequency
cat c | ngram 3 | cat> trigram.txt &
# Word frequency
cat d | ranked_frequency |cat > words.txt &
# scatter the stdin to the outher processes// nested tee block
cat a | tr -cs a-zA-Z \\n |  tee b c d > /dev/null &

cat teea | tee a e > /dev/null &

# the the file byte count
cat teeb | wc -c | cat > teea_res &
# redir the file to these two proc
cat ${IN} | tee teeb teea > /dev/null  &
sleep 1


# get the results and print everything? I am missing something here
cat e | 
  	sed 's/./&\
/g' |
	# Print absolute
	ranked_frequency |
	awk 'BEGIN {
		cat teea_res | getline NCHARS
		OFMT = "%.2g%%"}
		{print $1, $2, $1 / NCHARS * 100}' > acharacter.txt &
    rm a b c e f g d  trigram.txt words.txt digram.txt teea teeb teea_res teeb_res

