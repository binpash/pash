#!/bin/sh
# Automatically generated file
# Source file example/map-hierarchy.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Hierarchy map
# DESCRIPTION
# Given two directory hierarchies A and B passed as input arguments
# (where these represent a project at different parts of its lifetime)
# copy the files of hierarchy A to a new directory, passed as a third
# argument, corresponding to the structure of directories in B.
# Demonstrates the use of <em>join</em> to gather results from streams
# in the <em>scatter</em> part,
# and the use of <em>cat</em> to order asynchronous results in the gather part.
#
#  Copyright 2012-2013 Diomidis Spinellis
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

if [ ! -d "$1" -o ! -d "$2" -o -z "$3" ]
then
	echo "Usage: $0 dir-1 dir-2 new-dir-name" 1>&2
	exit 1
fi

NEWDIR="$3"

export LC_ALL=C

line_signatures()
{
	find $1 -type f -name '*.[chly]' -print |
	# Split path name into directory and file
	sed 's|\(.*\)/\([^/]*\)|\1 \2|' |
	while read dir file
	do
		# Print "directory filename content" of lines with
		# at least one alphabetic character
		# The fields are separated by ^A and ^B
		sed -n "/[a-z]/s|^|$dir$file|p" "$dir/$file"
	done |
	sort -t -k 2
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
	 {  line_signatures $1
}</dev/null  >$SGDIR/npfo-_unnamed_0.0
 {  line_signatures $2
}</dev/null  >$SGDIR/npfo-_unnamed_1.0
 {  join -t -1 2 -2 2 $SGDIR/npfo-_unnamed_0.0 $SGDIR/npfo-_unnamed_1.0 |
	# Print filename dir1 dir2
	sed 's///g' |
	awk -F 'BEGIN{OFS=" "}{print $1, $3, $4}' |
	# Unique occurrences
	sort -u
}</dev/null  >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
 {  awk '{print "mkdir -p \"'$NEWDIR'/" $3 "\""}' | sort -u
} <$SGDIR/npi-1.0.0 >$SGDIR/npfo-_unnamed_2.0
 {  awk '{print "cp \"" $2 "/" $1 "\" \"'$NEWDIR'/" $3 "/" $1 "\""}'
} <$SGDIR/npi-1.1.0 >$SGDIR/npfo-_unnamed_3.0

# Gather the results
	# Order: first make directories, then copy files
	./sgsh-tee -f -I -i $SGDIR/npfo-_unnamed_2.0 -i $SGDIR/npfo-_unnamed_3.0 |
	sh

)  3<&0 
