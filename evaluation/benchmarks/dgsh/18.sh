#!/bin/sh
# Automatically generated file
# Source file example/dir.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Directory listing
# DESCRIPTION
# Windows-like DIR command for the current directory.
# Nothing that couldn't be done with <code>ls -l | awk</code>.
# Demonstrates combined use of stores and streams.
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

FREE=`df -h . | awk '!/Use%/{print $4}'`

ls -n |
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
 {  awk '!/^total/ {print $6, $7, $8, $1, sprintf("%8d", $5), $9}'
} <$SGDIR/npi-0.0.0 >$SGDIR/npfo-_unnamed_0.0
NFILES=$(  {  wc -l
} <$SGDIR/npi-0.1.0 )
NDIRS=$(  {  grep -c '^d'
} <$SGDIR/npi-0.2.0 )
NBYTES=$(  {  awk '{s += $5} END {print s}'
} <$SGDIR/npi-0.3.0 )

# Gather the results
	cat $SGDIR/npfo-_unnamed_0.0
	echo "               `echo ${NFILES}` File(s) `echo ${NBYTES}` bytes"
	echo "               `echo ${NDIRS}` Dir(s) $FREE bytes free"

)  3<&0 
