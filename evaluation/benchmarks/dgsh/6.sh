#!/bin/sh
# Automatically generated file
# Source file example/word-properties.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Word properties
# DESCRIPTION
# Read text from the standard input and list words
# containing a two-letter palindrome, words containing
# four consonants, and words longer than 12 characters.
#
# Demonstrates the use of paste as a gather function
#
# Example:
# curl ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt | word-properties
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

IN=${VOC:-/usr/share/dict/words}
cat ${IN} | 
# Split input one word per line
tr -cs a-zA-Z \\n |
# Create list of unique words
sort -u |
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
	cat <&3 3<&-  >$SGDIR/npi-0.0.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.1.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.2.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.3.0
ln -s $SGDIR/npi-0.0.0 $SGDIR/npfo-_unnamed_0.0
 {  sed 's/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t
		g'
} <$SGDIR/npi-0.1.0 >$SGDIR/npfo-_unnamed_1.0
 {  sed -E 's/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t
		g'
} <$SGDIR/npi-0.2.0 >$SGDIR/npfo-_unnamed_2.0
 {  awk '{if (length($1) > 12) print "l:", length($1);
		else print ""}'
} <$SGDIR/npi-0.3.0 >$SGDIR/npfo-_unnamed_3.0

# Gather the results
	# Paste the four streams side-by-side
	paste $SGDIR/npfo-_unnamed_0.0 $SGDIR/npfo-_unnamed_1.0 $SGDIR/npfo-_unnamed_2.0 $SGDIR/npfo-_unnamed_3.0 |
	# List only words satisfying one or more properties
	grep :

)  3<&0 
