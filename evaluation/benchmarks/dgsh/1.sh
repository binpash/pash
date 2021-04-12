#!/bin/sh
# Automatically generated file
# Source file example/compress-compare.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS Compression benchmark
# DESCRIPTION
# Report file type, length, and compression performance for
# data received from the standard input.  The data never touches the
# disk.
# Demonstrates the use of stores.
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

IN=${FULL:-$PASH_TOP/evaluation/benchmarks/dgsh/input/dblp.xml}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input/}


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
	cat ${IN}  >$SGDIR/npi-0.0.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.1.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.2.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.3.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.4.0
NBYTES=$(  {  wc -c
} <$SGDIR/npi-0.0.0 )
FILETYPE=$(  {  file -
} <$SGDIR/npi-0.1.0 )
XZ=$(  {  xz -c | wc -c
} <$SGDIR/npi-0.2.0 )
BZIP2=$(  {  bzip2 -c | wc -c
} <$SGDIR/npi-0.3.0 )
GZIP=$(  {  gzip -c | wc -c
} <$SGDIR/npi-0.4.0 )

# Gather the results
	cat <<EOF
File type:	`echo ${FILETYPE}`
Original size:	`echo ${NBYTES}` bytes
gzip:		`echo ${GZIP}` bytes
bzip2:		`echo ${BZIP2}` bytes
xz:		`echo ${XZ}` bytes
EOF

)  3<&0 
