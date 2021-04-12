#!/bin/sh
# Automatically generated file
# Source file example/duplicate-files.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Find duplicate files
# DESCRIPTION
# List the names of duplicate files in the specified directory.
# Demonstrates the combination of streams with a relational join.
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

export IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input}

# Create list of files
find ${IN} -type f |

# Produce lines of the form
# MD5(filename)= 811bfd4b5974f39e986ddc037e1899e7
xargs openssl md5 |

# Convert each line into a "filename md5sum" pair
sed 's/^MD5(//;s/)= / /' |

# Sort by MD5 sum
sort -k2 |

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
 {  awk '{print $2}' | uniq -d
} <$SGDIR/npi-0.0.0 >$SGDIR/npfo-_unnamed_0.0
ln -s $SGDIR/npi-0.1.0 $SGDIR/npfo-_unnamed_1.0

# Gather the results
	# Join the repeated MD5 sums with the corresponding file names
	join -2 2 $SGDIR/npfo-_unnamed_0.0 $SGDIR/npfo-_unnamed_1.0 |
	# Output same files on a single line
	awk '
	BEGIN {ORS=""}
	$1 != prev && prev {print "\n"}
	END {if (prev) print "\n"}
	{if (prev) print " "; prev = $1; print $2}'

)  3<&0 
