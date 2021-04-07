#!/bin/bash
# tag: text properties
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#text-properties

set -e

IN=${VOC:-/usr/share/dict/words}
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

