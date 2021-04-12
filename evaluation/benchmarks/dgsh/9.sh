#!/bin/sh
# Automatically generated file
# Source file example/static-functions.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS C/C++ symbols that should be static
# DESCRIPTION
# Given as an argument a directory containing object files, show which
# symbols are declared with global visibility, but should have been
# declared with file-local (static) visibility instead.
# Demonstrates the use of streams and comm (1) to combine data from
# two sources.
#
#  Copyright 2014 Diomidis Spinellis
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


IN=${BIN:-/usr/local/bin}
# Find object files
find ${IN} -type f -exec grep -IL . "{}" \; |

  
# Print defined symbols
xargs nm |

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
 {  awk 'NF == 3 && $2 ~ /[A-Z]/ {print $3}' | sort
} <$SGDIR/npi-0.0.0 >$SGDIR/npfo-_unnamed_0.0
 {  awk '$1 == "U" {print $2}' | sort
} <$SGDIR/npi-0.1.0 >$SGDIR/npfo-_unnamed_1.0

# Gather the results
	# Print exports that are not imported
	comm -23 $SGDIR/npfo-_unnamed_0.0 $SGDIR/npfo-_unnamed_1.0

)  3<&0 
