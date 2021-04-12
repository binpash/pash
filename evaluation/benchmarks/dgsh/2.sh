#!/bin/sh
# Automatically generated file
# Source file example/commit-stats.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Git commit statistics
# DESCRIPTION
# Process the git history, and list the authors and days of the week
# ordered by the number of their commits.
# Demonstrates streams and piping through a function.
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

# Works without change
# Order by frequency
forder()
{
	sort |
	uniq -c |
	sort -rn
}

git log --format="%an:%ad" --date=default "$@" |
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
 {  awk -F: '{print $1}' | forder
} <$SGDIR/npi-0.0.0 >$SGDIR/npfo-_unnamed_0.0
 {  awk -F: '{print substr($2, 1, 3)}' | forder
} <$SGDIR/npi-0.1.0 >$SGDIR/npfo-_unnamed_1.0

# Gather the results
	echo "Authors ordered by number of commits"
	cat $SGDIR/npfo-_unnamed_0.0
	echo "Days ordered by number of commits"
	cat $SGDIR/npfo-_unnamed_1.0

)  3<&0 
