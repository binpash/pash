#!/bin/sh
# Automatically generated file
# Source file example/spell-highlight.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Highlight misspelled words
# DESCRIPTION
# Highlight the words that are misspelled in the command's standard
# input.
# Demonstrates stream buffering, the avoidance of pass-through
# constructs to avoid deadlocks, and the use of named streams.
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


export IN=${FULL:-$PASH_TOP/evaluation/benchmarks/dgsh/input/dblp.xml}

export LC_ALL=C

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
	cat ${IN}  >$SGDIR/npi-0.1.0
ln $SGDIR/npi-0.1.0 $SGDIR/npi-0.2.0
 {  sort /usr/share/dict/words
}</dev/null  >$SGDIR/npfo-dict.0
 {  tr -cs A-Za-z \\n |
	tr A-Z a-z |
	sort -u |
	comm -23 - $SGDIR/npfo-dict.0
} <$SGDIR/npi-0.1.0 >$SGDIR/npfo-errors.0
 {  cat
} <$SGDIR/npi-0.2.0 >$SGDIR/npfo-text.0

# Gather the results
	fgrep -f $SGDIR/npfo-errors.0 -i --color -w -C 2 $SGDIR/npfo-text.0

)  3<&0 
