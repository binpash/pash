#!/bin/sh
# Automatically generated file
# Source file example/code-metrics.sh
#!/usr/bin/env sgsh
#
# SYNOPSIS C code metrics
# DESCRIPTION
# Process a directory containing C source code, and produce a summary
# of various metrics.
# Demonstrates nesting, commands without input.
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
	 {  find "$@" \( -name \*.c -or -name \*.h \) -type f -print0
}</dev/null  >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
FNAMELEN=$(  {  tr \\0 \\n |
		# Remove path
		sed 's|^.*/||' |
		# Maintain average
		awk '{s += length($1); n++} END {print s / n}'
} <$SGDIR/npi-1.0.0 )
 {  xargs -0 cat
} <$SGDIR/npi-1.1.0 >$SGDIR/npi-2.0.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.1.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.2.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.3.0
ln $SGDIR/npi-2.0.0 $SGDIR/npi-2.4.0
 {  sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" |
			cpp -P 2>/dev/null
} <$SGDIR/npi-2.0.0 >$SGDIR/npi-3.0.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.1.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.2.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.3.0
ln $SGDIR/npi-3.0.0 $SGDIR/npi-3.4.0
NSTRUCT=$(  {   egrep -c 'struct[ 	]*{|struct[ 	]*[a-zA-Z_][a-zA-Z0-9_]*[       ]*{'
} <$SGDIR/npi-3.0.0 )
NTYPEDEF=$(  {  grep -cw typedef
} <$SGDIR/npi-3.1.0 )
NVOID=$(  {  grep -cw void
} <$SGDIR/npi-3.2.0 )
NGETS=$(  {  grep -cw gets
} <$SGDIR/npi-3.3.0 )
IDLEN=$(  {  tr -cs 'A-Za-z0-9_' '\n' |
				sort -u |
				awk '/^[A-Za-z]/ { len += length($1); n++ } END {print len / n}'
} <$SGDIR/npi-3.4.0 )
CHLINESCHAR=$(  {  wc -lc | awk '{OFS=":"; print $1, $2}'
} <$SGDIR/npi-2.1.0 )
NCCHAR=$(  {  sed 's/#/@/g' |
			cpp -traditional -P 2>/dev/null |
			wc -c |
			awk '{OFMT = "%.0f"; print $1/1000}'
} <$SGDIR/npi-2.2.0 )
NCOMMENT=$(  {  egrep -c '/\*|//'
} <$SGDIR/npi-2.3.0 )
NCOPYRIGHT=$(  {  grep -ci copyright
} <$SGDIR/npi-2.4.0 )
 {  find "$@" -name \*.c -type f -print0
}</dev/null  >$SGDIR/npi-4.0.0
ln $SGDIR/npi-4.0.0 $SGDIR/npi-4.1.0
 {  tr \\0 \\n
} <$SGDIR/npi-4.0.0 >$SGDIR/npi-5.0.0
ln $SGDIR/npi-5.0.0 $SGDIR/npi-5.1.0
NCFILE=$(  {  wc -l
} <$SGDIR/npi-5.0.0 )
NCDIR=$(  {  sed 's,/[^/]*$,,;s,^.*/,,' | sort -u | wc -l
} <$SGDIR/npi-5.1.0 )
 {  xargs -0 cat
} <$SGDIR/npi-4.1.0 >$SGDIR/npi-6.0.0
ln $SGDIR/npi-6.0.0 $SGDIR/npi-6.1.0
CLINESCHAR=$(  {  wc -lc | awk '{OFS=":"; print $1, $2}'
} <$SGDIR/npi-6.0.0 )
 {  sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" |
			cpp -P 2>/dev/null
} <$SGDIR/npi-6.1.0 >$SGDIR/npi-7.0.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.1.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.2.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.3.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.4.0
ln $SGDIR/npi-7.0.0 $SGDIR/npi-7.5.0
NFUNCTION=$(  {  grep -c '^{'
} <$SGDIR/npi-7.0.0 )
NGOTO=$(  {  grep -cw goto
} <$SGDIR/npi-7.1.0 )
NREGISTER=$(  {  grep -cw register
} <$SGDIR/npi-7.2.0 )
NMACRO=$(  {  grep -c '@[ 	]*define[ 	][ 	]*[a-zA-Z_][a-zA-Z0-9_]*('
} <$SGDIR/npi-7.3.0 )
NINCLUDE=$(  {  grep -c '@[ 	]*include'
} <$SGDIR/npi-7.4.0 )
NCONST=$(  {  grep -o -h '[0-9][x0-9][0-9a-f]*' | wc -l
} <$SGDIR/npi-7.5.0 )
NHFILE=$(  {  find "$@" -name \*.h -type f | wc -l
}</dev/null  )

# Gather the results
cat <<EOF
FNAMELEN: `echo ${FNAMELEN}`
NSTRUCT: `echo ${NSTRUCT}`
NTYPEDEF: `echo ${NTYPEDEF}`
IDLEN: `echo ${IDLEN}`
CHLINESCHAR: `echo ${CHLINESCHAR}`
NCCHAR: `echo ${NCCHAR}`
NCOMMENT: `echo ${NCOMMENT}`
NCOPYRIGHT: `echo ${NCOPYRIGHT}`
NCFILE: `echo ${NCFILE}`
NCDIR: `echo ${NCDIR}`
CLINESCHAR: `echo ${CLINESCHAR}`
NFUNCTION: `echo ${NFUNCTION}`
NGOTO: `echo ${NGOTO}`
NREGISTER: `echo ${NREGISTER}`
NMACRO: `echo ${NMACRO}`
NINCLUDE: `echo ${NINCLUDE}`
NCONST: `echo ${NCONST}`
NVOID: `echo ${NVOID}`
NHFILE: `echo ${NHFILE}`
NGETS: `echo ${NGETS}`
EOF

)  3<&0 
