#!/bin/sh
# Automatically generated file
# Source file example/web-log-stats.sh
#!/usr/bin/env sgsh -s /bin/bash
#
# SYNOPSIS Web log statistics
# DESCRIPTION
# Provides continuous statistics over web log stream data.
# Demonstrates stream processing.
# Provide as an argument either the name of a growing web log file
# or -s and a static web log file, which will be processed at a rate
# of about 10 lines per second.
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

# Size of the window to report in seconds
WINDOW=10
WINDOW_OLD=$(expr $WINDOW \* 2)

# Update interval in seconds
UPDATE=2

# Print the sum of the numbers read from the standard input
sum()
{
	awk '{ sum += $1 } END {print sum}'
}

# Print the rate of change as a percentage
# between the first (old) and second (new) value
change()
{
	# Can't use bc, because we have numbers in scientific notation
	awk "END {OFMT=\"%.2f%%\"; print ($2 - $1) * 100 / $1}" </dev/null
}

if [ "$1" = "-s" ]
then
	# Simulate log lines coming from a file
	while read line
	do
		echo $line
		sleep 0.01
	done  <"$2"
else
	tail -f "$1"
fi |
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
page=$(  {  awk -Winteractive '{print $7}'
} <$SGDIR/npi-0.0.0 )
 {  awk -Winteractive '{print $10}'
} <$SGDIR/npi-0.1.0 >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.2.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.3.0
total_bytes=$(  {  awk -Winteractive '{ s += $1; print s}'
} <$SGDIR/npi-1.0.0 )
total_pages=$(  {  awk -Winteractive '{print ++n}'
} <$SGDIR/npi-1.1.0 )
bytes=$(  { 
} <$SGDIR/npi-1.2.0 )
bytes_old=$(  { 
} <$SGDIR/npi-1.3.0 )

# Gather the results
	# Produce periodic reports
	while :
	do
		WINDOW_PAGES=$(echo ${bytes} -c | wc -l)
		WINDOW_BYTES=$(echo ${bytes} -c | sum )
		WINDOW_PAGES_OLD=$(echo ${bytes_old} -c | wc -l)
		WINDOW_BYTES_OLD=$(echo ${bytes_old} -c | sum)
		clear
		cat <<EOF
Total
-----
Pages: $(echo ${total_pages} -c)
Bytes: $(echo ${total_bytes} -c)

Over last ${WINDOW}s
--------------------
Pages: $WINDOW_PAGES
Bytes: $WINDOW_BYTES
kBytes/s: $(awk "END {OFMT=\"%.0f\"; print $WINDOW_BYTES / $WINDOW / 1000}" </dev/null )
Top page: $(echo ${page} -c | sort | uniq -c | sort -rn | head -1)

Change
------
Requests: $(change $WINDOW_PAGES_OLD $WINDOW_PAGES)
Data bytes: $(change $WINDOW_BYTES_OLD $WINDOW_BYTES)
EOF
	sleep $UPDATE
	done

)  3<&0 
