#!/bin/sh
# Automatically generated file
# Source file example/text-properties.sh
#!/usr/bin/env sgsh
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

# Consitent sorting across machines
export LC_ALL=C
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
IN=${VOC:-/usr/share/dict/words}
# Convert input into a ranked frequency list
ranked_frequency()
{
	awk '{count[$1]++} END {for (i in count) print count[i], i}' |
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

(

	export SGDIR=/tmp/sg-$$.0

	rm -rf $SGDIR

	# Cleanup on exit or interrupt
	cleanup()
	{
		SIGNAL=$1
		[ $SIGNAL = EXIT ] || echo sgsh interrupted. Cleaning up... 1>&2

		# Stop key-value stores
		
		# Kill processes we have launched in the background
		kill $SGPID 2>/dev/null

		# Remove temporary directory
		rm -rf "$SGDIR"

		# Propagate real signals and exit with non-0
		if [ $SIGNAL != EXIT ]
		then
			trap - $SIGNAL EXIT
			kill -s $SIGNAL $$
		fi

		# Exit with the original exit value
		exit

	}

	for sig in HUP INT QUIT TERM EXIT
	do
		trap "cleanup $sig" $sig
	done

	mkdir $SGDIR
	cat < ${IN}  >$SGDIR/npi-0.0.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.1.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.2.0
 {  tr -cs a-zA-Z \\n
} <$SGDIR/npi-0.0.0 >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.2.0
 {  ngram 2 >${OUT}/digram.txt
: ; } <$SGDIR/npi-1.0.0 >$SGDIR/npfo-none-1.0.0
 {  ngram 3 >${OUT}/trigram.txt
: ; } <$SGDIR/npi-1.1.0 >$SGDIR/npfo-none-1.1.0
 {  ranked_frequency >${OUT}/words.txt
: ; } <$SGDIR/npi-1.2.0 >$SGDIR/npfo-none-1.2.0
NCHARS=$(  {  wc -c
} <$SGDIR/npi-0.1.0 )
 {  sed 's/./&\
/g' |
	   # Print absolute and percentage
	   ranked_frequency |
	   awk 'BEGIN {OFMT = "%.2g%%"}
	   {print $1, $2, $1 / '"`echo ${NCHARS}`"' * 100}' >${OUT}/character.txt
: ; } <$SGDIR/npi-0.2.0 >$SGDIR/npfo-none-0.2.0

# Gather the results
./sgsh-tee  -i $SGDIR/npfo-none-1.0.0  -i $SGDIR/npfo-none-1.1.0  -i $SGDIR/npfo-none-1.2.0  -i $SGDIR/npfo-none-0.2.0 >/dev/null

)  3<&0 
