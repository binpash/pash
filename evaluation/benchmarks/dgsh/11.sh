#!/bin/sh
# FIXME
# Automatically generated file
# Source file example/committer-plot.sh
#!/usr/bin/env sgsh -s /bin/bash
#
# SYNOPSIS Plot git committer activity over time
# DESCRIPTION
# Process the git history, and create two PNG diagrams depicting
# committer activity over time. The most active committers appear
# at the center vertical of the diagram.
# Demonstrates image processing and no-output scatter blocks.
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


# Commit history in the form of ascending Unix timestamps, emails
git log --pretty=tformat:'%at %ae' |
awk 'NF == 2 && $1 > 100000 && $1 < '`date +%s` |
sort -n |
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
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.4.0
ln $SGDIR/npi-0.0.0 $SGDIR/npi-0.5.0
NCOMMITTERS=$(  {  awk '{print $2}' | sort -u | wc -l
} <$SGDIR/npi-0.0.0 )
LAST=$(  {  tail -1 | awk '{print $1}'
} <$SGDIR/npi-0.1.0 )
FIRST=$(  {  head -1 | awk '{print $1}'
} <$SGDIR/npi-0.2.0 )
NDAYS=$(  {  expr \( `echo ${LAST}` - `echo ${FIRST}` \) / 60 / 60  / 24
}</dev/null  )
 {  awk '{print $2}' |
	   sort |
	   uniq -c |
	   sort -n |
	   awk 'BEGIN {l = 0; r = '`echo ${NCOMMITTERS}`';}
		      {print NR % 2 ? l++ : --r, $2}' |
	   sort -k2
} <$SGDIR/npi-0.4.0 >$SGDIR/npfo-_unnamed_0.0
 {  sort -k2 |
	   join -j 2 - $SGDIR/npfo-_unnamed_0.0 |
	   # Order by time
	   sort -k 2n |
	{
	# Create portable bitmap
	echo 'P1'
	echo "`echo ${NCOMMITTERS}` `echo ${NDAYS}`"
	perl -na -e '
	BEGIN { @empty['`echo ${NCOMMITTERS}`' - 1] = 0; @committers = @empty; }
	sub out { print join("", map($_ ? "1" : "0", @committers)), "\n"; }

	$day = int($F[1] / 60 / 60 / 24);
	$pday = $day if (!defined($pday));

	while ($day != $pday) {
		out();
		@committers = @empty;
		$pday++;
	}

	$committers[$F[2]] = 1;

	END { out(); }
	'
	} |
	# Enlarge points into discs through morphological convolution
	pgmmorphconv -erode <(
cat <<EOF
P1
7 7
0 0 0 1 0 0 0
0 0 1 1 1 0 0
0 1 1 1 1 1 0
1 1 1 1 1 1 1
0 1 1 1 1 1 0
0 0 1 1 1 0 0
0 0 0 1 0 0 0
EOF
	)
} <$SGDIR/npi-0.5.0 >$SGDIR/npi-1.0.0
ln $SGDIR/npi-1.0.0 $SGDIR/npi-1.1.0
 {  pnmtopng >large.png
: ; } <$SGDIR/npi-1.0.0 >$SGDIR/npfo-none-1.0.0
 {  pamscale -width 640 |
		   pnmtopng >small.png
: ; } <$SGDIR/npi-1.1.0 >$SGDIR/npfo-none-1.1.0

# Gather the results
./sgsh-tee  -i $SGDIR/npfo-none-1.0.0  -i $SGDIR/npfo-none-1.1.0 >/dev/null

)  3<&0 
